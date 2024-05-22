cd("/Users/julianashwin/Documents/GitHub/EcoNNet.jl")

"""
Add extra cores if you want to parallelise later
"""

using LaTeXStrings, TableView, CSV
using Distributed
#addprocs(2)
#rmprocs(2)
workers()

@everywhere include("src/EcoNNet_dev.jl")

save_gifs = false

"""
Define some options.
Expectations are denoted by a capital E prefix,
leads and lags are specified as a `_lead` or `_lag` suffix
"""

@everywhere options = EcoNNetOptions(infoset = [:y_lag, :ϵ_π, :ϵ_y],
    expectations = [:Eπ,:Ey,:Eπ_lead],
    endogenous = [:π, :y],
    exogenous = [:ϵ_π,:ϵ_y],
    states = [:y_lag, :ϵ_π, :ϵ_y],
	auxiliary = [:r],
    N = 500000, num_nodes = 32, activation = relu, window = 40000,
	burnin = 50000, learning_gap = 10000, plotting_gap = 10000)

@everywhere beliefs = initialise_beliefs(options)

"""
Define the parameters as a Named Tuple.
"""

@everywhere par = (β = 0.95, κ = 0.05, η = 0.75, σ = 0.25,
	ϕ_π = 0.5, π_star = 1.0, α = 0.75,
	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);

"""
State the equilibrium conditions of the model as a function which returns
    a vector of zeros
"""

@everywhere function equilibrium_conditions_fast(F::Array{Float64,1},
    x::Array{Float64,1},states_input::Array{Float64,1},predictions_input::Array{Float64,1})
    # Manually unpack the states
    #p_lag::Float64 = states_input[1]
    y_lag::Float64 = states_input[1]
    ϵ_π::Float64 = states_input[2]
	ϵ_y::Float64 = states_input[3]
    # and the predictions
    #Ep::Float64 = predictions_input[1]
    #Ey::Float64 = predictions_input[2]
    Eπ_lead::Float64 = predictions_input[3]
    # and the endogenous variables
    π_t::Float64 = x[1]
    y_t::Float64 = x[2]

	# NKPC
    F[1] = par.β*Eπ_lead + par.κ*y_t + ϵ_π - π_t ::Float64

	if π_t > par.π_star
		r_t = par.ϕ_π*π_t + par.α*(π_t - par.π_star)
	elseif π_t < -par.π_star
		r_t = par.ϕ_π*π_t + par.α*(π_t + par.π_star)
	else
		r_t = par.ϕ_π*π_t
	end

	# IS curve
    F[2] = par.η*y_lag - par.σ*(r_t - Eπ_lead) + ϵ_y - y_t ::Float64

    return F
end


"""
Define equations of motion under perfect foresight as a useful comparison
"""
function perfect_foresight(inputs)
    # Manually unpack the states (endogenous variables in same order as in options)
    π_lag = inputs[1]
    y_lag = inputs[2]
    # and the predictions
    # π[t+1] = 1/β*(π[t] - κ*y[t])
    π = (1/par.β)*(π_lag - par.κ*y_lag)
    # y[t+1] = η*y[t] - σ*(ϕ_π*π[t+1] + α π[t+1]^3 - π[t+2])

	if π > par.π_star
		r = par.ϕ_π*π + par.α*(π - par.π_star)
	elseif π < -par.π_star
		r = par.ϕ_π*π + par.α*(π + par.π_star)
	else
		r = par.ϕ_π*π
	end

	y = (par.β/(par.β+par.σ*par.κ))*(
		par.η*y_lag - par.σ*r + (par.σ/par.β)*π)
	# Impose upper and lower bounds to allow plotting
	π = min(π,max(π,-1e6),1e6)
	y = min(y,max(y,-1e6),1e6)

    outputs = [π,y]

end


"""
Define the variables and expectations to keep track of.
All variables which appear as an expectation need to be included here.
Lagged values which appear as state variables need not be included.
"""

@everywhere variables = Symbol.(cat(Vector(options.exogenous),
    Vector(options.endogenous),outputnames(options),
    Vector(options.expectations),Vector(options.auxiliary), dims = 1));
s = DataFrame(ones(options.N, length(variables)), variables);
@everywhere indices = extract_indices(options, variables);


"""
Plot steady state conditions
"""
function NKPC_condition(y)
    π = par.κ/(1-par.β)*y
end
function Taylor_condition(π)
	if π > par.π_star
		r = par.ϕ_π*π + par.α*(π - par.π_star)
	elseif π < -par.π_star
		r = par.ϕ_π*π + par.α*(π + par.π_star)
	else
		r = par.ϕ_π*π
	end
end
function IS_condition(π)
    y = -(par.σ/(1-par.η))*(Taylor_condition(π) - π)
end

pyplot()
plot_points = -4.0:0.01:4.0;
function plot_ss(plot_points)
	# Set up plot
	ss_plot = plot(xlabel = L"y_{t-1}", xlims = (minimum(plot_points),maximum(plot_points)),
    	ylabel = L"\pi_t", ylims = (minimum(plot_points),maximum(plot_points)),legend=:bottomright, yguidefontrotation=-90)
	# Plot NKPC condition
	plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
	# Plot IS condition
	plot!(IS_condition.(plot_points),plot_points, label = "IS Curve", color = :green)
	# Fill in indeterminate area
	plot!(plot_points, par.π_star*ones(len((plot_points))),
		fillrange=[-par.π_star*ones(len(plot_points))], fillalpha = 0.5,
		 color = :paleturquoise1, label = L"\pi_t < \pi^*")
	display(plot!(plot_points, -par.π_star*ones(len((plot_points))),
		color = :paleturquoise1, label = false))
	return ss_plot
end

#par = @set par.π_star = 1.0
plot_ss(plot_points)

#par = @set par.η = 0.95
#par = @set par.σ = 0.25
#par = @set par.κ = 0.05
#plot_ss(plot_points)



"""
Find the steady states numerically
"""
# Define a steady state function to numerically solve
"""
Function that gives zero when x (inflation, output) is at a deterministic steady state
"""
function steady_states(F::Array{Float64,1},x::Array{Float64,1})
    π::Float64 = x[1]
    y::Float64 = x[2]

    F[1]::Float64 = par.β*π + par.κ*y - π
    F[2]::Float64 = par.η*y - par.σ*(Taylor_condition(π) - π) - y
    return F
end
# Exogenous processes
ss = Dict{Symbol,Float64}()
ss[:ϵ_π] = 0.0
ss[:ϵ_y] = 0.0
# upper steady state
sstates = nlsolve(steady_states, [2.0, 2.0])
upper = deepcopy(ss);
upper[:π] = sstates.zero[1];
upper[:y] = sstates.zero[2];
upper[:Eπ_lead] = upper[:π];
# central steady state
sstates = nlsolve(steady_states, [0.0, 0.0])
central = deepcopy(ss);
central[:π] = sstates.zero[1];
central[:y] = sstates.zero[2];
central[:Eπ_lead] = central[:π];
# lower steady state
sstates = nlsolve(steady_states, [-2.0, -2.0])
lower = deepcopy(ss);
lower[:π] = sstates.zero[1];
lower[:y] = sstates.zero[2];
lower[:Eπ_lead] = lower[:π];

"""
Plot perfect foresight paths
"""
#for tt in 1:50
initial_ss = deepcopy(central)
starts = [(π=-0.3,y=-3.5,periods=100,arrows=[10,25,98]),
	(π=0.3,y=3.5,periods=100,arrows=[10,25,98]),
	(π=-2.8,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=2.8,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-2.2,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=2.2,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-1.898,y=-3.5,periods=80,arrows=[10,55]),
	(π=-2.394,y=-3.5,periods=60,arrows=[25]),
	(π=1.898,y=3.5,periods=80,arrows=[10,55]),
	(π=2.394,y=3.5,periods=60,arrows=[25])
	]
starts = [(π=-0.61,y=-3.5,periods=30,arrows=[10,25]),
	(π=0.61,y=3.5,periods=30,arrows=[10,25]),
	(π=-0.,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=0.,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-0.9,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=0.9,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-2.0,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=2.0,y=3.5,periods=100,arrows=[10,50,98]),
	(π=1.5,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=-1.5,y=3.5,periods=100,arrows=[10,50,98])
	]
plot_ss(plot_points)
for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods])
	if start == starts[1]
		phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 1:(start[:periods]-1), label = "Perfect foresight paths", arrow_size = .5)
	else
		phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 1:(start[:periods]-1), label = "", arrow_size = .5)
	end
end
plot!(size = (600,400))

## Save these paths to a file
all_paths_df = DataFrame(zeros(1,10), 
	[:ϵ_π, :ϵ_y, :π, :y, :Eπ, :Ey, :Eπ_lead, :r, :y_lag, :path_number])
for ii in 1:length(starts)
	start = starts[ii];
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods]);
	paths1[:,:y_lag] = vcat([NaN],paths1.y[1:(start.periods-1)]);
	paths1[:,:path_number] .= ii
	all_paths_df = vcat(all_paths_df, paths1)
end
all_paths_df = all_paths_df[all_paths_df.path_number .!= 0.,:]
CSV.write("figures/pw_linear/sim_data/perf_foresight_paths_eta0p75.csv", all_paths_df)

	

#savefig("figures/pw_linear/perf_phase_pistar"*rep_pnt(par.π_star)*#
#	"_alpha"*rep_pnt(par.α)*".pdf")


"""
Test the equilibrium conditions and step! function by confirming the steady states
"""
# Steady state should give zeros from equilibrium_conditions_fast
s = initialise_df(s, upper);
t = 3;
if len(indices.outputindex_current)>0
	predictions1 = cat(Array(s[t-1, indices.outputindex_current]), Array(s[t-1,indices.outputindex_lead]), dims = 1);
else
	predictions1 = cat(Array(s[t-1,indices.outputindex_lead]), dims = 1);
end
states1 = extract_states(s,t,indices);
starting_values1 = Vector(s[t,indices.endogindex]);
F1 = zeros(length(options.endogenous))
display(step_fast!(cat(starting_values1, states1, predictions1, dims = 1),options))
display(equilibrium_conditions_fast(F1,starting_values1, states1, predictions1))

"""
Initialise beliefs by training on (some of the) steady state(s)
"""
# Initialise s
s = initialise_df(s, lower, gap = 500, steadystate_alt = upper)
s = initialise_df(s, central)
# Initialise and train beliefs
@everywhere beliefs = initialise_beliefs(options)
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss)


"""
Run learning simulation
"""
# Training options
options.burnin_use_net = true;
options.burnin = 100000;
options.learning_gap = 100000;
options.plotting_gap = 100000;
options.window = 99999;
options.plot_vars = [:π, :y, :Eπ, :Ey]

# Simulate the learning for a set number of periods
noise_π = par.σ_π*randn((options.N - options.burnin))
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:];
s.ϵ_π[(options.burnin+1):options.N] = simulate_ar(par.ρ_π, par.σ_π, options.N - options.burnin, noise_π)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs,s = simulate_learning(options.burnin:options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-5999):(500000-4999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.π[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))
#savefig("figures/pw_linear/sim_series_pistar"*rep_pnt(par.π_star)*
#	"_alpha"*rep_pnt(par.α)*".pdf")

export_df = deepcopy(s)
#export_df = export_df[options.N-99999:options.N,:]
rename!(export_df, replace.(names(export_df), "π" => "pi"))
rename!(export_df, replace.(names(export_df), "ϵ" => "epsilon"))
export_df.r = Taylor_condition.(export_df.pi)
#export_df = export_df[:,[:epsilon_pi, :epsilon_y, :pi, :y, :r, :y_lag, :Epi_lead]]

CSV.write("figures/pw_linear/sim_data/export_data_eta0p75.csv", export_df)

#save("networks/pwlin_eta"*rep_pnt(par.η)*".jld2", "beliefs", beliefs)
#beliefs = load("networks/pwlin_eta"*rep_pnt(par.η)*".jld2", "beliefs", "beliefs");

"""
Run a whole bunch of simulations
"""
## Initialise DF to add to
plm_df = DataFrame(zeros(1,7), 
	[:version, :progress, :y_lag, :y, :pi, :Epi_lead, :Rsq_Epi])
gr() # Set GR backend for plots as it's the fastest
for jj in 1:10
	print("Initialising version "*string(jj)*"\n")
	s = initialise_df(s, central)
	# Initialise and train beliefs
	@everywhere beliefs = initialise_beliefs(options)
	@time beliefs = learn!(beliefs, s, options.N, options, indices, loss)
	# Run learning updates a bunch of times
	for ii in 1:50
		# Keep last version of previous run as burnin
		s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:];
		# Random shock draws
		noise_π = par.σ_π*randn((options.N - options.burnin))
		s.ϵ_π[(options.burnin+1):options.N] = simulate_ar(par.ρ_π, par.σ_π, options.N - options.burnin, noise_π)
		noise_y = par.σ_y*randn((options.N - options.burnin))
		s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
		# Simulate distribution and update beliefs
		beliefs,s = simulate_learning(options.burnin:options.N, s, beliefs, indices, options)
		Rsq_π = 1 - var(s.π-s.Eπ)/var(s.π)
		# Get results
		y_opts = -3.5:0.01:3.5
		run_df = DataFrame(version = jj.*ones(length(y_opts)), 
							progress = ii.*ones(length(y_opts)), 
							y_lag = y_opts, 
							y = zeros(length(y_opts)),
							pi = zeros(length(y_opts)),
							Epi_lead = zeros(length(y_opts)), 
							Rsq_Epi = Rsq_π.*ones(length(y_opts)))
		for yy in 1:nrow(run_df)
			inputs = [run_df.y_lag[yy], 0., 0.]
			predictions = predict!(inputs, beliefs)
			run_df.pi[yy] = predictions[1]
			run_df.y[yy] = predictions[2]
			run_df.Epi_lead[yy] = predictions[3]
		end
		plm_df = vcat(plm_df, run_df)
	end
	plm_df = plm_df[plm_df.version .!= 0.,:]
	CSV.write("figures/pw_linear/sim_data/many_sims_eta0p"*string(Int(par.η*100))*".csv", plm_df)
end
		
groups = string.(plm_df.version).*"-".*string.(plm_df.progress)
plot(plm_df.y_lag, plm_df.Epi_lead, group = groups)




"""
Assess the accuracy of solution
"""
s.Eπ_error = s.Eπ-s.π
s.Ey_error = s.Ey-s.y
s[:,:y_lag] .= 0.0
s.y_lag[2:options.N] = s.y[1:(options.N-1)]
s[:,:π_lag] .= 0.0
s.π_lag[2:options.N] = s.π[1:(options.N-1)]

s.Eπ_lead_error = vcat((s.Eπ_lead[1:(options.N-1)]-s.π[2:(options.N)]),[0.])
plot(s.Eπ_lead_error[(options.N-499):options.N])

plot_range=410001:500000
plot(layout = (1,2))
histogram!(s.π[plot_range], label = "Inflation", subplot = 1, legend = :bottomright)
histogram!(s.y[plot_range], label = "Output", color = :red, subplot = 2, legend = :bottomright)

# Calculate the R squared
Rsq_π = 1 - var(s.π-s.Eπ)/var(s.π)
Rsq_y = 1 - var(s.y-s.Ey)/var(s.y)

# Conduct DHM test
pvals = DHM_test(s, 5000:500:500000, 400, hvars = [:y_lag], include_const = true);
pvals = DHM_test(s, 5000:500:500000, 400, hvars = [:π_lag], include_const = true);
pvals = DHM_test(s, 5000:500:500000, 400, hvars = [:ϵ_π], include_const = true);
pvals = DHM_test(s, 5000:500:500000, 400, hvars = [:ϵ_y], include_const = true);
pvals = DHM_test(s, 5000:500:500000, 400, hvars = [:y_lag, :ϵ_π, :ϵ_y], include_const = true);


"""
Create a heatmap
"""
y_range = Array(-6.0:0.05:6.0)
π_range= Array(-3.5:0.05:3.5)
heatmap_df = deepcopy(s[plot_range,:])
heatmap_df.y_lag = s[plot_range.-1,:y]
heatmap_df.y_grid = map(x -> findnearest(y_range,x), Array(heatmap_df.y_lag))
heatmap_df.π_grid = map(x -> findnearest(π_range,x), Array(heatmap_df.π))

heatmap_mat = zeros(len(π_range),len(y_range))
prog = Progress(len(π_range), dt = 1, desc="Creating heatmap: ")
for ππ in 1:len(π_range)
	for yy in 1:len(y_range)
		instances = (heatmap_df[:,:π_grid].==π_range[ππ]).* (heatmap_df[:,:y_grid].==y_range[yy])
		heatmap_mat[ππ,yy] = sum(instances)
	end
	next!(prog)
end

heatmap(y_range,π_range,heatmap_mat, c=:coolwarm,
	xlabel= L"y_{t-1}",ylabel=L"\pi_t")
plot!(size=(800,300))
#savefig("figures/heatmap_illus_sim.pdf")


avg_error = zeros(len(y_range))
for yy in 1:len(y_range)
	short_df = heatmap_df[1:89999,:]
	instances = Array(short_df[:,:y_grid].==y_range[yy])
	errors = short_df[instances,:Eπ_lead_error]
	avg_error[yy] = mean(errors)
	if isnan(avg_error[yy])
		avg_error[yy] = 0.0
	end
end
# Plot the distribution of forecast errors
pyplot()
plot(layout = (2,1), link = :x)
heatmap!(y_range,π_range,heatmap_mat, subplot=1, c=[:white,:blue],
	ylabel=L"\pi_{t}",colorbar =:none,yguidefontrotation=-90)
plot!(y_range, avg_error,subplot=2, xlabel= L"y_{t-1}", ylabel= L"Avg. [E \pi_{t+1} - \pi_{t+1}]",
	legend=:false,yguidefontrotation=0)
plot!(size=(600,400))
savefig("figures/heatmap_errors_pwlin_sim.pdf")


using HypothesisTests

rrr = 490000:500000
CorrelationTest(s.Eπ_lead_error[rrr],s.y[rrr.-1])
CorrelationTest(s.Eπ_lead_error[rrr],s.ϵ_π[rrr])
CorrelationTest(s.Eπ_lead_error[rrr],s.ϵ_y[rrr])
CorrelationTest(s.Eπ_lead_error[rrr.-1],s.ϵ_π[rrr])
pvalue(CorrelationTest(s.Eπ_lead_error[rrr],s.y[rrr.-1]))





"""
Identify stochastic steady states
"""
# Upper steady state
upper_stoch = deepcopy(upper)
paths = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 1000,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false)
upper_stoch[:π] = paths.π[1000];
upper_stoch[:Eπ_lead] = paths.Eπ_lead[1000];
upper_stoch[:y] = paths.y[1000]

# Lower steady state
lower_stoch = deepcopy(lower)
paths = irf(:ϵ_y, lower_stoch, beliefs, shock_period = 5, periods = 1000,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false)
lower_stoch[:π] = paths.π[1000];
lower_stoch[:Eπ_lead] = paths.Eπ_lead[1000];
lower_stoch[:y] = paths.y[1000]

pyplot()
paths1 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 100,
	magnitude = -0.0, persistence = par.ρ_y, show_plot = true,plot_vars=[:y,:π])
paths2 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 100,
	magnitude = -0.03, persistence = par.ρ_y, show_plot = true,
	plot_vars=[:y,:π,:Ey, :Eπ])






"""
Plot phase diagram
"""

pyplot()
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=3.0,y=3.0,periods=100,arrows=[10]),
	(π=-3.,y=-3.,periods=100,arrows=[10]),
	(π=0.5,y=0.5,periods=100,arrows=[22, 65]),
	(π=-0.5,y=-0.5,periods=100,arrows=[22, 65]),
	(π=-.2,y=-0.2,periods=100,arrows=[19]),
	(π=0.3,y=0.2,periods=100,arrows=[19]),
	#(π=1.0,y=-1.0,periods=100,arrows=[9]),
	#(π=-1.0,y=-1.0,periods=100,arrows=[9])
	]
for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = start[:periods],
		magnitude = 0.0, persistence = par.ρ_y, show_plot = false)
	paths.y_lag = cat(start[:y],paths.y[1:(start.periods -1 )],dims=1)
	if start == starts[1]
		phase_arrow_plot(paths, [:y_lag,:π], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 2:(start[:periods]), label = "Equilibrium paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths, [:y_lag,:π], arrow_points=start[:arrows].-5, h_points = 2:(start[:periods]),
			v_points = 2:start[:periods], arrow_size = .5, final_arrow = true)
		end
end
plot!(size = (600,400))
#savefig("figures/pw_linear/phase_nnet_pistar"*rep_pnt(par.π_star)*
#	"_alpha"*rep_pnt(par.α)*".pdf")


for yy in -3.5:0.5:3.5
	inputs = [yy, 0., 0.]
	predictions = predict!(inputs, beliefs)[3]
	print("y[t-1] = "*string(yy)*", Eπ[t+1] = "*string(round(predictions, digits = 4))*"\n")
end

for yy in -3.5:0.5:3.5
	inputs = [yy, 0., 0.]
	predictions = predict!(inputs, beliefs)[1]
	print("y[t-1] = "*string(yy)*", Eπ[t] = "*string(round(predictions, digits = 4))*"\n")
end