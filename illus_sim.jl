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


"""
Define some options.
Expectations are denoted by a capital E prefix,
leads and lags are specified as a `_lead` or `_lag` suffix
"""

@everywhere options = EcoNNetOptions(infoset = [:π_lag, :y_lag, :ϵ_π, :ϵ_y],
    expectations = [:Eπ,:Ey,:Eπ_lead],
    endogenous = [:π, :y],
    exogenous = [:ϵ_π,:ϵ_y],
    states = [:y_lag, :ϵ_π, :ϵ_y],
	auxiliary = [:r],
    N = 500000, num_nodes = 24, activation = σ, window = 40000)

@everywhere beliefs = initialise_beliefs(options)

"""
Define the parameters as a Named Tuple.
"""

@everywhere par = (β = 0.95, κ = 0.05,
	η = 0.95, σ = 0.25, ϕ_π = 0.5,
    α_2 = 0.0, α_3 = 0.075,
	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);

@everywhere par = (β = 0.95, κ = 0.05,
	η = 0.95, σ = 0.25, ϕ_π = 1.5,
    α_2 = 0.0, α_3 = 0.0,
	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);


@everywhere par = (β = 1.051, κ = -0.0258,
	η = 0.8505, σ = 0.3879, ϕ_π = 0.5326,
    α_2 = 0.0, α_3 = 0.0,
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
    π::Float64 = x[1]
    y::Float64 = x[2]

    # p[t] = ϕ_pp*pe[t+1] + ϕ_py*y[t] + α_2*y[t]^2 - α_3*y[t]^3 + ϵ_p[t]
    F[1] =  par.β*Eπ_lead + par.κ*y + ϵ_π - π ::Float64
    # y[t] = η*y[t-1] + ϕ_yp*p[t] + ϵ_y[t]
	r = par.ϕ_π*π + par.α_3*π^3
    F[2] = par.η*y_lag - par.σ*(r - Eπ_lead) + ϵ_y - y ::Float64

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
    πt = (1/par.β)*(π_lag - par.κ*y_lag)
    # y[t+1] = η*y[t] - σ*(ϕ_π*π[t+1] + α π[t+1]^3 - π[t+2])
	y = (par.β/(par.β+par.σ*par.κ))*(
		par.η*y_lag - par.σ*(((par.ϕ_π*par.β-1)/par.β)*πt + par.α_3*πt^3))
	# Impose upper and lower bounds to allow plotting
	πt = min(πt,max(πt,-1e100),1e100)
	y = min(y,max(y,-1e100),1e100)

    outputs = [πt,y]

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
Calculate steady states
"""
function NKPC_condition(y)
    π = par.κ/(1-par.β)*y
end
function IS_condition(π)
    y = (par.σ/(1-par.η))*((1-par.ϕ_π)*π - par.α_3*π^3)
end
function steady_states(F::Array{Float64,1},x::Array{Float64,1})
    π::Float64 = x[1]
    y::Float64 = x[2]
    F[1]::Float64 = par.β*π + par.κ*y - π
    F[2]::Float64 = par.η*y - par.σ*(par.ϕ_π*π + par.α_3*π^3 - π) - y
    return F
end

# Exogenous processes
ss = Dict{Symbol,Float64}()
ss[:ϵ_π] = 0.0
ss[:ϵ_y] = 0.0
# upper steady state
sstates = nlsolve(steady_states, [2.0, 2.0])
upper = Dict{Symbol,Float64}();
upper[:π] = sstates.zero[1];
upper[:y] = sstates.zero[2];
upper[:Eπ_lead] = upper[:π];
# central steady state
sstates = nlsolve(steady_states, [0.0, 0.0])
central = Dict{Symbol,Float64}();
central[:π] = sstates.zero[1];
central[:y] = sstates.zero[2];
central[:Eπ_lead] = central[:π];
# lower steady state
lower = Dict{Symbol,Float64}();
sstates = nlsolve(steady_states, [-2.0, -2.0])
lower[:π] = sstates.zero[1];
lower[:y] = sstates.zero[2];
lower[:Eπ_lead] = lower[:π];

s = initialise_df(s, upper);
s = initialise_df(s, ss);
s = initialise_df(s, lower, gap = 2, steadystate_alt = upper)


"""
Plot steady state conditions and perfect foresight paths
"""

@everywhere par = (β = 1.06, κ = -0.03,
	η = 0.85, σ = 0.41, ϕ_π = 0.55,
    α_2 = 0.0, α_3 = 0.0,
	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);

@everywhere par = (β = 0.0882, κ = 0.6971,
	η = 0.9390, σ = 0.1335, ϕ_π = 0.5754,
    α_2 = 0.0, α_3 = 0.0,
	ρ_y = 0.4933, σ_y = 0.1983, ρ_π = 0.5981, σ_π = 0.3806);

#using Plotly, PlotlyJS
pyplot()
plot_points = -4.0:0.01:4.0;
function plot_ss(plot_points)
	ss_plot = plot(xlabel = L"y_{t-1}", xlims = (-4.0,4.0),
    	ylabel = L"\pi_t", ylims = (-4.0, 4.0),legend=:bottomright, yguidefontrotation=-90)
		plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
		display(plot!(IS_condition.(plot_points),plot_points, label = "IS Curve", color = :green))
		# Plot perfect foresight paths
	return ss_plot
end

#for tt in 1:50
starts = [(π=0.75,y=-1.5,periods=100,arrows=[10,60]),
	(π=-0.75,y=1.5,periods=100,arrows=[10,60]),
	(π=0.0,y=-1.5,periods=100,arrows=[10,60]),
	(π=0.0,y=1.5,periods=100,arrows=[10,60]),
	(π=-1.024,y=-3.0,periods=80,arrows=[10,55]),
	(π=-2.350785,y=-3.0,periods=60,arrows=[25]),
	(π=1.024,y=3.0,periods=80,arrows=[10,55]),
	(π=2.350785,y=3.0,periods=60,arrows=[25])]
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [
(π=2.2,y=3.0,periods=3,arrows=[2]),
(π=0.74,y=1.0,periods=3,arrows=[2]),
(π=2.3,y=3.0,periods=3,arrows=[2]),
(π=0.8,y=1.0,periods=3,arrows=[2]),
(π=-2.2,y=-3.0,periods=3,arrows=[2]),
(π=-0.74,y=-1.0,periods=3,arrows=[2]),
(π=-2.3,y=-3.0,periods=3,arrows=[2]),
(π=-0.8,y=-1.0,periods=3,arrows=[2])]
for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods], lags = 1)

	if start == starts[1]
		phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
			v_points = 1:(start[:periods]), label = "Perfect foresight paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
			v_points = 1:(start[:periods]), label = "", arrow_size = .5, final_arrow = true)
	end
end
plot!(size = (600,400))
plot!(title = "Perfect foresight estimated Determinate model")
savefig("figures/phase_illus_perf.pdf")

if false
	starts = [(π=0.756541,y=1.0,periods=10,arrows=[1])]

	periods = 5

	anim = @animate for tt in 2:periods
		plot_ss(plot_points)
		plot!(size = (600,400))
		#plot!(xlims = (-30,30))
		for start in starts
			initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
			paths1 = pf_path(initial_ss, periods = periods)


			if start == starts[1]
				phase_arrow_plot(paths1, [:y,:π], arrow_points=[], h_points = 2:min(start[:periods], tt),
					v_points = 1:(min(start[:periods], tt)-1), label = "Perfect foresight paths", arrow_size = .5)
				if tt >=3
					phase_arrow_plot(paths1, [:y,:π], arrow_points=[min(start[:periods], tt)-2], h_points = 2:min(start[:periods], tt),
						v_points = 1:(min(start[:periods], tt)-1), label = "", arrow_size = .5)
				end
			else
				phase_arrow_plot(paths1, [:y,:π], arrow_points=[], h_points = 2:min(start[:periods], tt),
					v_points = 1:(min(start[:periods], tt)-1), label = "", arrow_size = .5)
				if tt >=3
					phase_arrow_plot(paths1, [:y,:π], arrow_points=[min(start[:periods], tt)-2], h_points = 2:min(start[:periods], tt),
						v_points = 1:(min(start[:periods], tt)-1), label = "", arrow_size = .5)
				end
			end
		end

	end
	gif(anim, "figures/phase_illus_perf.gif", fps = 2)
end


"""
Test the equilibrium conditions and step! function by confirming the steady states
"""

s = initialise_df(s, lower);
s = initialise_df(s, ss);
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

@everywhere beliefs = initialise_beliefs(options)
s = initialise_df(s, lower, gap = 500, steadystate_alt = upper)
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss)


"""
Run learning simulation
"""
noise_π = par.σ_π*randn(nrow(s))
s.ϵ_π = simulate_ar(par.ρ_π, par.σ_π, options.N, noise_π)
noise_y = par.σ_y*randn(nrow(s))
s.ϵ_y = simulate_ar(par.ρ_y, par.σ_y, options.N, noise_y)
plot(s.ϵ_y[1:200])
options.burnin = 50000;
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], lower, gap = 500, steadystate_alt = upper);
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
options.burnin_use_net = false;
options.window = 49999;
options.learning_gap = 50000;
options.plotting_gap = 50000;
options.plot_vars = [:π, :y, :Eπ, :Ey]

# Simulate the learning for a set number of periods
noise_π = par.σ_π*randn((options.N - options.burnin))
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_π[(options.burnin+1):options.N] = simulate_ar(par.ρ_π, par.σ_π, options.N - options.burnin, noise_π)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs,s = simulate_learning((options.burnin+1):options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-4999):(500000-0)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.π[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))
savefig("figures/illus_sim_series.pdf")

export_df = s[options.N-99999:options.N,:]
rename!(export_df, replace.(names(export_df), "π" => "pi"))
rename!(export_df, replace.(names(export_df), "ϵ" => "epsilon"))
export_df.r = par.ϕ_π.*export_df.pi .+ par.α_3.*export_df.pi.^3

CSV.write("estimation/illustrative/illus_sim.csv", export_df)


"""
Plot phase diagram
"""

pyplot()
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=3.0,y=3.0,periods=100,arrows=[10]),
	(π=-3.,y=-3.,periods=100,arrows=[10]),
	(π=-0.008,y=-0.008,periods=100,arrows=[42, 65]),
	(π=-0.05,y=-0.05,periods=100,arrows=[42, 65]),
	#(π=1.0,y=1.0,periods=100,arrows=[9]),
	#(π=-1.0,y=1.0,periods=100,arrows=[9]),
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
savefig("figures/phase_illus_sim.pdf")

starts = [(π=3.0,y=3.0,periods=100,arrows=[10]),
	(π=-3.,y=-3.,periods=100,arrows=[10]),
	(π=-0.008,y=-0.008,periods=100,arrows=[40, 65]),
	(π=-0.03,y=-0.03,periods=100,arrows=[40, 65]),
	#(π=1.0,y=1.0,periods=100,arrows=[9]),
	#(π=-1.0,y=1.0,periods=100,arrows=[9]),
	#(π=1.0,y=-1.0,periods=100,arrows=[9]),
	#(π=-1.0,y=-1.0,periods=100,arrows=[9])
	]
periods = 100
anim = @animate for tt in 2:periods
	plot_ss(plot_points)
	for start in starts
		initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
		paths1 = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = tt,
			magnitude = 0.0, persistence = par.ρ_y, show_plot = false)
		paths1.y_lag = cat(start[:y],paths1.y[1:(tt -1 )],dims=1)
		if start == starts[1]
			phase_arrow_plot(paths1, [:y_lag,:π], arrow_points=[], h_points = 2:min(start[:periods], tt),
				v_points = 2:(min(start[:periods], tt)), label = "Perfect foresight paths", arrow_size = .5)
			if tt >=3
				phase_arrow_plot(paths1, [:y_lag,:π], arrow_points=[min(start[:periods], tt)-2], h_points = 2:min(start[:periods], tt),
					v_points = 2:(min(start[:periods], tt)), label = "", arrow_size = .5)
			end
		else
			phase_arrow_plot(paths1, [:y_lag,:π], arrow_points=[], h_points = 2:min(start[:periods], tt),
				v_points = 2:(min(start[:periods], tt)), label = "", arrow_size = .5)
			if tt >=3
				phase_arrow_plot(paths1, [:y_lag,:π], arrow_points=[min(start[:periods], tt)-2], h_points = 2:min(start[:periods], tt),
					v_points = 2:(min(start[:periods], tt)), label = "", arrow_size = .5)
			end
		end
	end

end
gif(anim, "figures/phase_illus_sim.gif", fps = 15)


"""
Identify stochastic steady states
"""
# Upper steady state
upper_stoch = deepcopy(upper)
paths = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 1000,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false)
upper_stoch[:π] = paths.π[1000];
upper_stoch[:Eπ_lead] = paths.Eπ_lead[1000];
upper_stoch[:y] = paths.y[1000];

# Lower steady state
lower_stoch = deepcopy(lower)
paths = irf(:ϵ_y, lower_stoch, beliefs, shock_period = 5, periods = 1000,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false)
lower_stoch[:π] = paths.π[1000];
lower_stoch[:Eπ_lead] = paths.Eπ_lead[1000];
lower_stoch[:y] = paths.y[1000];


paths1 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 50,
	magnitude = -0.5, persistence = par.ρ_y, show_plot = false)
paths2 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 50,
	magnitude = 0.5, persistence = par.ρ_y, show_plot = false)
paths3 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 50,
	magnitude = -1.5, persistence = par.ρ_y, show_plot = false)
paths4 = irf(:ϵ_y, lower_stoch, beliefs, shock_period = 5, periods = 50,
	magnitude = 1.5, persistence = par.ρ_y, show_plot = false)
plot(paths1.π, label =L"\pi_t", xlabel = "Periods", legend = :bottomright,ylims = (-3.0,3.0))
plot!(paths1.y, label =L"y_t");plot!(paths1.ϵ_y, label =L"\epsilon_{y,t}")
plot!(size=(300,200))
savefig("figures/irf1_illus_sim.pdf")
plot(paths2.π, label =L"\pi_t", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths2.y, label =L"y_t");plot!(paths2.ϵ_y, label =L"\epsilon_{y,t}")
plot!(size=(300,200))
savefig("figures/irf2_illus_sim.pdf")
plot(paths3.π, label =L"\pi_t", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths3.y, label =L"y_t");plot!(paths3.ϵ_y, label =L"\epsilon_{y,t}")
plot!(size=(300,200))
savefig("figures/irf3_illus_sim.pdf")
plot(paths4.π, label =L"\pi_t", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths4.y, label =L"y_t");plot!(paths4.ϵ_y, label =L"\epsilon_{y,t}")
plot!(size=(300,200))
savefig("figures/irf4_illus_sim.pdf")

# All plots together
irf_plot = plot(layout = (2,2),ylims = (-3.0,3.0),size=(1000,600))
plot!(paths1.π, label =L"\pi_t", xlabel = "Periods", legend = :bottomright, subplot=1)
plot!(paths1.y, label =L"y_t");plot!(paths1.ϵ_y, label =L"\epsilon_{y,t}")
plot!(paths2.π, label =L"\pi_t", xlabel = "Periods", legend = false, subplot=2)
plot!(paths2.y, label =L"y_t", subplot=2);plot!(paths2.ϵ_y, label =L"\epsilon_{y,t}", subplot=2)
plot!(paths3.π, label =L"\pi_t", xlabel = "Periods", legend = false, subplot=3)
plot!(paths3.y, label =L"y_t", subplot=3);plot!(paths3.ϵ_y, label =L"\epsilon_{y,t}", subplot=3)
plot!(paths4.π, label =L"\pi_t", xlabel = "Periods", legend = false, subplot=4)
plot!(paths4.y, label =L"y_t", subplot=4);plot!(paths4.ϵ_y, label =L"\epsilon_{y,t}", subplot=4)
plot!(size=(600,400))
display(irf_plot)
savefig("figures/irfs_illus_sim.pdf")



"""
Compare PLM and ALM
"""
s.Eπ_error = s.Eπ-s.π
s.Ey_error = s.Ey-s.y
s.Eπ_lead_error = vcat((s.Eπ_lead[1:(options.N-1)]-s.π[2:(options.N)]),[0.])
plot(s.Eπ_lead_error)

plot_range=410001:500000
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
savefig("figures/heatmap_illus_sim.pdf")


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


pyplot()
plot(layout = (2,1), link = :x)
heatmap!(y_range,π_range,heatmap_mat, subplot=1, c=[:white,:blue],
	ylabel=L"\pi_{t}",colorbar =:none,yguidefontrotation=-90)
plot!(y_range, avg_error,subplot=2, xlabel= L"y_{t-1}", ylabel= L"Avg. [E \pi_{t+1} - \pi_{t+1}]",
	legend=:false,yguidefontrotation=0)
plot!(size=(600,400))
savefig("figures/heatmap_errors_illus_sim.pdf")



Rsq_π = 1 - var(heatmap_df.π-heatmap_df.Eπ)/var(heatmap_df.π)
Rsq_y = 1 - var(heatmap_df.y-heatmap_df.Ey)/var(heatmap_df.y)


"""
End of script
"""
