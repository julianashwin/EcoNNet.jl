cd("/Users/julianashwin/Documents/DPhil/nnet_learning_julia")

"""
Add extra cores if you want to parallelise later
"""
using Distributed
addprocs(4)
workers()

@everywhere include("EcoNNet/src/EcoNNet.jl")


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
    N = 500000, num_nodes = 64, activation = σ, window = 40000)

@everywhere beliefs = initialise_beliefs(options)

"""
Define the parameters as a Named Tuple.
"""

@everywhere par = (β = 0.95, κ = 0.05,
	η = 0.95, σ = 0.25, ϕ_π = 0.5,
    α_2 = 0.0, α_3 = 0.075,
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
    π = (1/par.β)*(π_lag - par.κ*y_lag)
    # y[t+1] = η*y[t] - σ*(ϕ_π*π[t+1] + α π[t+1]^3 - π[t+2])
	y = (par.β/(par.β+par.σ*par.κ))*(
		par.η*y_lag - par.σ*(((par.ϕ_π*par.β-1)/par.β)*π + par.α_3*π^3))
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
#using Plotly, PlotlyJS
#plotly()
plot_points = -4.0:0.01:4.0;
ss_plot = plot(xlabel = "Output (state)", xlims = (-4.0,4.0),
    ylabel = "Inflation (control)", ylims = (-4.0, 4.0),legend=:topright)
plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
display(plot!(IS_condition.(plot_points),plot_points, label = "IS Curve", color = :green))
# Plot perfect foresight paths
initial_ss = deepcopy(central)
starts = [(π=0.75,y=-1.5,periods=100,arrows=[10,60]),
	(π=-0.75,y=1.5,periods=100,arrows=[10,60]),
	(π=0.0,y=-1.5,periods=61,arrows=[10,60]),
	(π=0.0,y=1.5,periods=61,arrows=[10,60]),
	(π=-1.024,y=-3.0,periods=100,arrows=[10,55]),
	(π=-2.351,y=-3.0,periods=50,arrows=[25]),
	(π=1.024,y=3.0,periods=100,arrows=[10,55]),
	(π=2.351,y=3.0,periods=50,arrows=[25])]
for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods])
	phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
		v_points = 1:start[:periods])
end
display(ss_plot)
savefig("figures/phase_illus_perf.png")

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

s = initialise_df(s, lower, gap = 500, steadystate_alt = upper)
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss)


"""
Run learning simulation
"""
noise_y = par.σ_y*randn(nrow(s))
s.ϵ_y = simulate_ar(par.ρ_y, par.σ_y, options.N, noise_y)
plot(s.ϵ_y[1:2000])
options.burnin = 50000;
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], lower, gap = 500, steadystate_alt = upper);
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
options.burnin_use_net = false;
options.learning_gap = 500;
options.plotting_gap = 500;
options.plot_vars = [:π, :y, :Eπ, :Ey]

# Simulate the learning for a set number of periods
gr() # Set GR backend for plots as it's the fastest
s[300000:400000,:]= s[400000:options.N,:]
@time beliefs,s = simulate_learning(100000:options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (400000-3999):(400000-2999)
plot(layout=(2,1), xlabel = "Periods",legend = false)
plot!(s.π[plot_range], subplot = 1, ylabel = "Inflation")
plot!(s.y[plot_range], subplot = 2, ylabel = "Output")
plot!(size = (1000,400), left_margin = 7mm)
savefig("figures/illus_sim_series.png")


"""
Plot phase diagram
"""

phase_plot = plot(xlabel = "Output (state)", xlims = (-4.0,4.0),
	ylabel = "Inflation (control)", ylims = (-4.0, 4.0),legend = :topright)
plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
plot!(IS_condition.(plot_points),plot_points, label = "IS Curve", color = :green)

initial_ss = deepcopy(central)
starts = [(π=3.0,y=3.0,periods=100,arrows=[10,80]),
	(π=-3.,y=-3.,periods=100,arrows=[10,80]),
	(π=-0.1,y=-0.1,periods=110,arrows=[30,40,100]),
	(π=-0.05,y=-0.08,periods=140,arrows=[45,55,135]),
	#(π=1.0,y=1.0,periods=100,arrows=[9]),
	#(π=1.0,y=1.0,periods=100,arrows=[9]),
	#(π=-1.0,y=-1.0,periods=100,arrows=[9]),
	#(π=-1.0,y=-1.0,periods=100,arrows=[9])
	]
for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths = irf(:ϵ_y, initial_ss, beliefs, shock_period = 5, periods = start[:periods],
		magnitude = 0.0, persistence = par.ρ_y, show_plot = false)
	paths.y_lag = cat(start[:y],paths.y[1:(start.periods -1 )],dims=1)
	phase_arrow_plot(paths, [:y_lag,:π], arrow_points=start[:arrows].-5, h_points = 5:(start[:periods]-1),
		v_points = 6:start[:periods])
end
display(phase_plot)
savefig("figures/phase_illus_sim.png")


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
plot(paths1.π, label ="Inflation", xlabel = "Periods", legend = :bottomright,ylims = (-3.0,3.0))
plot!(paths1.y, label ="Output");plot!(paths1.ϵ_y, label ="Shock")
savefig("figures/irf1_illus_sim.png")
plot(paths2.π, label ="Inflation", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths2.y, label ="Output");plot!(paths2.ϵ_y, label ="Shock")
savefig("figures/irf2_illus_sim.png")
plot(paths3.π, label ="Inflation", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths3.y, label ="Output");plot!(paths3.ϵ_y, label ="Shock")
savefig("figures/irf3_illus_sim.png")
plot(paths4.π, label ="Inflation", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths4.y, label ="Output");plot!(paths4.ϵ_y, label ="Shock")
savefig("figures/irf4_illus_sim.png")

# All plots together
irf_plot = plot(layout = (2,2),ylims = (-3.0,3.0),size=(1000,600))
plot!(paths1.π, label ="Inflation", xlabel = "Periods", legend = :topright, subplot=1)
plot!(paths1.y, label ="Output");plot!(paths1.ϵ_y, label ="Shock")
plot!(paths2.π, label ="Inflation", xlabel = "Periods", legend = false, subplot=2)
plot!(paths2.y, label ="Output", subplot=2);plot!(paths2.ϵ_y, label ="Shock", subplot=2)
plot!(paths3.π, label ="Inflation", xlabel = "Periods", legend = false, subplot=3)
plot!(paths3.y, label ="Output", subplot=3);plot!(paths3.ϵ_y, label ="Shock", subplot=3)
plot!(paths4.π, label ="Inflation", xlabel = "Periods", legend = false, subplot=4)
plot!(paths4.y, label ="Output", subplot=4);plot!(paths4.ϵ_y, label ="Shock", subplot=4)

display(irf_plot)
savefig("figures/irf_illus_sim.png")



"""
Compare PLM and ALM
"""
s.Eπ_error = s.Eπ-s.π
s.Ey_error = s.Ey-s.y
s.Eπ_lead_error = vcat((s.Eπ_lead[1:(options.N-1)]-s.π[2:(options.N)]),[0.])
plot(s.Eπ_lead_error)

plot_range=310001:500000
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
heatmap(y_range,π_range,heatmap_mat,c=ColorGradient([:white,:blue]),
	xlabel= "Output (state)",ylabel="Inflation (control)")
plot!(size=(800,300))
savefig("figures/heatmap_illus_sim.png")


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


plot(layout = (2,1), link = :x)
heatmap!(y_range,π_range,heatmap_mat, subplot=1, c=ColorGradient([:white,:blue]),
	ylabel="Inflation (control)",legend =:left)
plot!(y_range, avg_error,subplot=2, xlabel= "Output (state)", ylabel= "Avg. Inflation forecast error",
	left_margin = 5mm,legend=:false)
savefig("figures/heatmap_errors_illus_sim.png")



Rsq_π = 1 - var(heatmap_df.π-heatmap_df.Eπ)/var(heatmap_df.π)
Rsq_y = 1 - var(heatmap_df.y-heatmap_df.Ey)/var(heatmap_df.y)


"""
Use Network in grid solution method
"""
p_y = (1 + par.ρ_y)/2
q_y = p_y
N = 101
ψ_y = sqrt(par.σ_y^2/(1 - par.ρ_y^2))*sqrt(N-1)
ϵ_y_range, ϵ_y_kernel = Rouwenhorst(p_y,q_y,ψ_y,N)
gr()
heatmap(ϵ_y_range,
    ϵ_y_range, ϵ_y_kernel,
    c=cgrad([:white, :blue]),
    xlabel="shock(t+1)", ylabel="shock(t)",
    title="My title")


display(join(["Order in grid should be ",indices.statenames_all]))
#p_range = -5.0:0.25:5.0
y_range = -3.5:0.1:3.5
#ϵ_y_range = Array{Float64}(-2.5:0.25:2.5)  # par.σ_y* randn(10)
#ϵ_y_density = round.(pdf.(Normal(0.0, par.σ_y),ϵ_y_range), digits = 12)
#ϵ_y_density = len(ϵ_y_density).*ϵ_y_density./sum(ϵ_y_density)

# Create state_grid
state_grid = (vec(collect(Iterators.product(y_range, ϵ_y_range))))
state_grid = [(collect(x)) for x in state_grid]
state_grid = permutedims(reshape(hcat(state_grid...), (length(state_grid[1]), length(state_grid))))
#@eval @everywhere state_grid=$state_grid

println(join(["Order in shock range should be:",(indices.statenames_current)]," "))
shock_range = (vec(collect(Iterators.product(ϵ_y_range))))
shock_range = [(collect(x)) for x in shock_range]
shock_range = permutedims(reshape(hcat(shock_range...), (length(shock_range[1]), length(shock_range))))

ϵ_y_transition = DataFrame(ϵ_y = repeat(ϵ_y_range, inner = length(ϵ_y_range)),
	ϵ_y_lead = repeat(ϵ_y_range, outer = length(ϵ_y_range)))
ϵ_y_transition.prob = vec(transpose(ϵ_y_kernel))

# Create a tuple with the transition probabilities (if only one shock remember to include comma)
transition_probabilities = (ϵ_y = ϵ_y_transition,)



# Create_df_grid makes DataFrame with states and expectations for each grid point
@time df_grid = create_df_grid(state_grid, beliefs, indices)
# step_1_map fills in endogenous variables for period t
df_grid = step_1_map(df_grid, beliefs, indices, options, sval = 0.0)
# Create the transition probability for each point
@time df_grid_new, loss_weights, probs_df =
	create_trans_probs(df_grid, transition_probabilities, indices,
		shock_range, parallel = true)

# step_2_prep fills in predictions, shocks and their probability for t+1
@time grid_map, fill_cols = step_2_prep(df_grid_new, beliefs, indices, options, sval = 0.0)
# in_grid_map will be used to solve for t+1 realisations in parallel, so needs to be available to each worker
in_grid_map  = SharedArray{Float64}(grid_map)
# out_grid_map will contain the t+1 realisations for each state and t+1 shock
out_grid_map  = SharedArray{Float64}(size(grid_map,1), len(indices.endognames))
# Solve forward using distributed for (first time will compile so be a lot slower so maybe pre-run on subset)
@time begin
	@sync @distributed for ii in 1:size(in_grid_map,1)
		out_grid_map[ii,:] = step_fast!(in_grid_map[ii,:],options)
	end
end

# df_grid_new can now be used to train the network
df_grid_new[:,fill_cols] = out_grid_map

beliefs = update_beliefs(df_grid_new, beliefs, indices, options,
	epochs = 2, cycles = 100, verbose = true, weights = loss_weights)

"""
Iteratively update beliefs based on the state-space grid
"""
gr()
for it in 1:10
    display(println("Running through iteration ", it))
	global df_grid = create_df_grid(state_grid, beliefs, indices)
	df_grid = step_1_map(df_grid, beliefs, indices, options, sval = 0.0)
	@time df_grid_new, loss_weights, probs_df = create_trans_probs(df_grid, transition_probabilities,
		indices, shock_range, parallel = true)
	# step_2_prep fills in predictions, shocks and their probability for t+1
	@time grid_map, fill_cols = step_2_prep(df_grid_new, beliefs, indices, options, sval = 0.0)
	# Grid_map will be used to solve for t+1 realisations in parallel, so needs to be available to each worker
	# in_grid_map will be used to solve for t+1 realisations in parallel, so needs to be available to each worker
	global in_grid_map  = SharedArray{Float64,2}(grid_map)
	# out_grid_map will contain the t+1 realisations for each state and t+1 shock
	global out_grid_map = SharedArray{Float64,2}(size(grid_map,1), len(indices.endognames))

	@time begin
		@sync @distributed for ii in 1:size(in_grid_map,1)
			out_grid_map[ii,:] = step_fast!(in_grid_map[ii,:],options)
		end
	end
	# df_grid_new can now be used to train the network
	df_grid_new[:,fill_cols] = out_grid_map
	@time global beliefs = update_beliefs(df_grid_new, beliefs, indices, options,
		epochs = 1, cycles = 200, verbose = true, weights = loss_weights)

	starting_point = deepcopy(upper)
	starting_point[:ϵ_y] = 0.0
	paths, plt1 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	    magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	    plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.5,3.5])
	starting_point = deepcopy(lower)
	starting_point[:ϵ_y] = 0.0
	paths, plt2 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	    magnitude = 0.0, persistence = par.ρ_y, show_plot = true,
	    plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.5,3.5])
	display(plot(plt1, plt2, layout = (2,1), title = join(["Iteration: ", it])))
end



display(join(["Order in grid should be ",indices.statenames_all]))
y_range = Array(y_range)


ranges = (y_lag = y_range, ϵ_y_range = ϵ_y_range)
state_t_vars = [:y_lag, :ϵ_y]
state_tp1_vars = [:y, :ϵ_y_lead]


for it in 1:30
    display(println("Running through iteration ", it))
	global df_grid = create_df_grid(state_grid, beliefs, indices)
	df_grid.state_t = 1:nrow(df_grid)
	df_grid.state_tp1 = Int.(zeros(nrow(df_grid)))
	df_grid.uncon_prob = zeros(nrow(df_grid))
	@time df_grid = step_1_map(df_grid, beliefs, indices, options, sval = 0.0)

	# Make sure that the step1 variables correspond to grid points
	df_grid.y = map(x -> findnearest(ranges.y_lag,x), Array(df_grid.y))

	# Find step 2 points and transition probabilities
	@time global df_grid_new, loss_weights, probs_df = create_trans_probs(df_grid, transition_probabilities,
		indices, shock_range, parallel = true)
	# step_2_prep fills in predictions, shocks and their probability for t+1
	@time grid_map, fill_cols = step_2_prep(df_grid_new, beliefs, indices, options, sval = 0.0)
	# Grid_map will be used to solve for t+1 realisations in parallel, so needs to be available to each worker
	# in_grid_map will be used to solve for t+1 realisations in parallel, so needs to be available to each worker
	global in_grid_map  = SharedArray{Float64,2}(grid_map)
	# out_grid_map will contain the t+1 realisations for each state and t+1 shock
	global out_grid_map = SharedArray{Float64,2}(size(grid_map,1), len(indices.endognames))

	@time begin
		@sync @distributed for ii in 1:size(in_grid_map,1)
			out_grid_map[ii,:] = step_fast!(in_grid_map[ii,:],options)
		end
	end
	# df_grid_new can now be used to train the network
	df_grid_new[:,fill_cols] = out_grid_map

	# Make sure that the step 2 variables correspond to grid points
	df_grid_new.y_lead = map(x -> findnearest(ranges.y_lag,x), Array(df_grid_new.y_lead))

	# Classify state for t+1
	poss_states = nrow(df_grid)
	prog = Progress(poss_states, dt = 1, desc="Classifying t+1 state: ")
	for st in 1:poss_states
		#state_defs = unique(df_grid_new[:,state_t_vars])
		state_def = df_grid[st,state_t_vars]
		sel_rows = BitArray(((df_grid_new.y .== state_def.y_lag).*
			(df_grid_new.ϵ_y_lead .== state_def.ϵ_y)))
		df_grid_new[sel_rows,:state_tp1] .= st
		next!(prog)
	end

	# Generate transition probability matrix for whole system
	@time global CP = create_kernel(df_grid_new, poss_states)
	#display(heatmap(CP))
	display("Finding eigenvectors")
	CP ./= sum(CP,dims = 2)

	global A_I = (transpose(CP)-I(size(CP,1)))
	#@elapsed MP_null = nullspace(A_I,atol=1e-15,rtol=1e-15)
	mc = MarkovChain(CP)
	@time MP = stationary_distributions(mc)[1]
	plot(1:length(MP),MP)

	MP = MP./sum(MP)
	#plot(1:len(MP),MP)

	# Use Bayes Rule to find unconditional probabilities
	MP_mat = repeat(MP, outer = [1, len(MP)])
	UP = MP_mat.*CP
	UP ./= sum(UP)

	# These unconditional probabilities need to be
	for st in 1:nrow(df_grid_new)
		df_grid_new[st,:uncon_prob] = UP[df_grid_new.state_t[st],df_grid_new.state_tp1[st]]
	end

	display(plot(df_grid_new.y, df_grid_new.uncon_prob))
	df_train = df_grid_new[(df_grid_new.uncon_prob .> 1e-16),:]

	loss_weights = Matrix(transpose(df_train.uncon_prob))
	loss_weights ./=(sum(loss_weights)/size(loss_weights,2))
	loss_weights = repeat(loss_weights, len(indices.expectnames_all))

	@time global beliefs = update_beliefs(df_train, beliefs, indices, options,
		epochs = 4, cycles = 500, verbose = true, weights = loss_weights)
	if true
		starting_point = deepcopy(upper)
		starting_point[:ϵ_y] = 0.0
		paths, plt1 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	    	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	    	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
		starting_point = deepcopy(lower)
		starting_point[:ϵ_y] = 0.0
		paths, plt2 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	    	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	    	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
		display(plot(plt1, plt2, layout = (2,1), title = join(["Iteration: ", it])))
	end
end
save("networks/meandyn_2eq_24nodes.jld2", "beliefs", beliefs)
beliefs = load("networks/meandyn_2eq_both_24nodes.jld2", "beliefs");

plot(df_grid_new.y, df_grid_new.uncon_prob, ylabel = "Probability", xlabel = "y")
savefig("figures/uncondist_2eq.png")
starting_point = deepcopy(upper)
starting_point[:ϵ_y] = 0.0
paths1, plt1 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
starting_point = deepcopy(lower)
starting_point[:ϵ_y] = 0.0
paths2, plt2 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
display(plot(plt1, plt2, layout = (2,1)))
savefig("figures/irf_2eq.png")

starting_point = deepcopy(upper)
starting_point[:y] = 3.5
paths1, plt1 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
starting_point = deepcopy(lower)
starting_point[:y] = -3.5
paths2, plt2 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
starting_point = deepcopy(central)
starting_point[:y] = 0.1
paths3, plt1 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
starting_point = deepcopy(central)
starting_point[:y] = -0.1
paths4, plt2 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])




phase_plot = plot(xlabel = "y",
    xlims = (-4.0,4.0),
    ylabel = "p",
    ylims = (-9.0, 9.0))
output = -4.0:0.01:4.0;
plot!(output,Δp_condition.(output), label = "p condition", color = "black")
plot!(output,Δy_condition.(output), label = "y condition", color = "black")

phase_arrow_plot(paths1, [:y,:p], arrow_points=[20,60], h_points = 11:99,
	v_points = 12:100)
phase_arrow_plot(paths2, vars, arrow_points=[2,20], h_points = 11:99,
	v_points = 12:100)
phase_arrow_plot(paths3, vars, arrow_points=[10], h_points = 11:99,
	v_points = 12:100)
phase_arrow_plot(paths4, vars, arrow_points=[20], h_points = 11:99,
	v_points = 12:100)
savefig("figures/phase_2eq_both.png")
