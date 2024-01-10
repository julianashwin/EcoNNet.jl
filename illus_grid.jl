cd("/Users/julianashwin/Documents/GitHub/EcoNNet.jl")

using Distributed
addprocs(2)
#rmprocs(2)
@everywhere include("src/EcoNNet_dev.jl")

"""
Define some options.
Expectations are denoted by a capital E prefix,
leads and lags are specified as a `_lead` or `_lag` suffix
"""

@everywhere options = EcoNNetOptions(infoset = [:y_lag, :ϵ_y],
    expectations = [:Eπ_lead],
    endogenous = [:π, :y],
    exogenous = [:ϵ_y],
    states = [:y_lag, :ϵ_y],
	auxiliary = [:r],
    N = 500000, num_nodes = 24, activation = σ, window = 40000)

@everywhere beliefs = initialise_beliefs(options)

"""
Define the parameters as a Named Tuple.
"""
# Parameters need to pe avaialble to all cores
@everywhere par = (β = 0.95, κ = 0.05,
	ϕ_yy = 0.95, σ = 0.25, ϕ_π = 0.5,
    α_2 = 0.0, α_3 = 0.075,
	ρ_y = 0.5, σ_y = 0.5, ρ_π=0.5, σ_π=0.3);

"""
State the equilibrium conditions of the model as a function which returns
    a vector of zeros
"""
# Equilibrium conditions need to avaialble to all cores
@everywhere function equilibrium_conditions_fast(F::Array{Float64,1},
    x::Array{Float64,1},states_input::Array{Float64,1},predictions_input::Array{Float64,1})
    # Manually unpack the states
    #π_lag::Float64 = states_input[1]
    y_lag::Float64 = states_input[1]
    ϵ_y::Float64 = states_input[2]
    # and the predictions
    #Ep::Float64 = predictions_input[1]
    #Ey::Float64 = predictions_input[2]
    Eπ_lead::Float64 = predictions_input[1]
    # and the endogenous variables
    π::Float64 = x[1]
    y::Float64 = x[2]

    # p[t] = ϕ_pp*pe[t+1] + ϕ_py*y[t] + α_2*y[t]^2 - α_3*y[t]^3 + ϵ_p[t]
    F[1] =  par.β*Eπ_lead + par.κ*y - π ::Float64
    # y[t] = ϕ_yy*y[t-1] + ϕ_yp*p[t] + ϵ_y[t]
	r = par.ϕ_π*π + par.α_3*π^3
    F[2] = par.ϕ_yy*y_lag - par.σ*(r - Eπ_lead) + ϵ_y - y ::Float64

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
    # y[t+1] = ϕ_yy*y[t] - σ*(ϕ_π*π[t+1] + α π[t+1]^3 - π[t+2])
	y = (par.β/(par.β+par.σ*par.κ*par.β))*(
		par.ϕ_yy*y_lag - par.σ*(((par.ϕ_π*par.β-1)/par.β)*π + par.α_3*π^3))
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
    y = (par.σ/(1-par.ϕ_yy))*((1-par.ϕ_π)*π - par.α_3*π^3)
end
function steady_states(F::Array{Float64,1},x::Array{Float64,1})
    π::Float64 = x[1]
    y::Float64 = x[2]
    F[1]::Float64 = par.β*π + par.κ*y - π
    F[2]::Float64 = par.ϕ_yy*y - par.σ*(par.ϕ_π*π + par.α_3*π^3 - π) - y
    return F
end

# Exogenous processes
ss = Dict{Symbol,Float64}()
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
gr()
plot_points = -4.0:0.01:4.0;
ss_plot = plot(xlabel = "Output (state)", xlims = (-4.0,4.0),
    ylabel = "Inflation (control)", ylims = (-4.0, 4.0))
plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
display(plot!(IS_condition.(plot_points),plot_points, label = "IS Curve", color = :green))

# Plot perfect foresight paths
initial_ss = deepcopy(central)
starts = [(π=0.75,y=-1.5,periods=100,arrows=[10,60]),
	(π=-0.75,y=1.5,periods=100,arrows=[10,60]),
	(π=0.0,y=-1.5,periods=61,arrows=[10,60]),
	(π=0.0,y=1.5,periods=61,arrows=[10,60]),
	(π=-1.05,y=-3.0,periods=100,arrows=[10,55]),
	(π=-2.356,y=-3.0,periods=50,arrows=[25]),
	(π=1.05,y=3.0,periods=100,arrows=[10,55]),
	(π=2.356,y=3.0,periods=50,arrows=[25])]
for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods])
	phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
		v_points = 1:start[:periods])
end
#plot!(size = (1000,600))
display(ss_plot)
#savefig("figures/phase_illus_perf.png")
gr()
"""
Test the equilibrium conditions and step! function by confirming the steady states
"""

s = initialise_df(s, upper);
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
s = initialise_df(s, upper, gap = 500, steadystate_alt = lower)
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss,
	n_cycles=100)


"""
Run learning simulation
"""
noise_y = par.σ_y*randn(nrow(s))
s.ϵ_y = simulate_ar(par.ρ_y, par.σ_y, options.N, noise_y)
plot(s.ϵ_y[1:2000])
options.burnin = 50000;
#s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], lower, gap = 500, steadystate_alt = upper);
#s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
options.burnin_use_net = false;
options.learning_gap = 100;
options.plotting_gap = 100;
options.window = 20000;
options.plot_vars = [:π, :y, :Eπ_lead]

# Simulate the learning for a set number of periods
@time beliefs,s = simulate_learning(options.burnin:(options.burnin+100000), s, beliefs,
	indices, options)
tt = options.burnin + 100000
plot(legend = :bottomleft, xguidefontsize=8)
plot!(s[(tt - options.plotting_window):tt,:π], label = "Inflation")
plot!(s[(tt - options.plotting_window):tt,:y], label = "Output")
plot!(s[(tt - options.plotting_window):tt,:Eπ_lead], label = "Expected Inflation")
display(plot!(title = ""))
#savefig("figures/realtime_learning.png")




"""
Use Network in grid solution method
"""
p_y = (1 + par.ρ_y)/2
q_y = p_y
N = 21
ψ_y = sqrt(par.σ_y^2/(1 - par.ρ_y^2))*sqrt(N-1)
ϵ_y_range, ϵ_y_kernel = Rouwenhorst(p_y,q_y,ψ_y,N)
ϵ_y_range, ϵ_y_kernel = Rouwenhorst(p_y,q_y,N)
gr()
heatmap(ϵ_y_range,
    ϵ_y_range, ϵ_y_kernel,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="shock(t+1)", ylabel="shock(t)",
    title="My title")


display(join(["Order in grid should be ",indices.statenames_all]))
π_range = Array(-4.0:0.1:4.0)
y_range = Array(-6.0:0.1:6.0)
#ϵ_y_range = Array{Float64}(-2.5:0.25:2.5)  # par.σ_y* randn(10)
#ϵ_y_density = round.(pdf.(Normal(0.0, par.σ_y),ϵ_y_range), digits = 12)
#ϵ_y_density = len(ϵ_y_density).*ϵ_y_density./sum(ϵ_y_density)

# Create state_grid
state_labels = (vec(collect(Iterators.product(y_range, ϵ_y_range))))
state_labels = [(collect(x)) for x in state_labels]
state_grid = permutedims(reshape(hcat(state_labels...), (length(state_labels[1]), length(state_labels))))
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

"""
Mean Dynamics
"""

π_dict = Dict(π_range[i] => i for i = 1:len(π_range))
y_dict = Dict(y_range[i] => i for i = 1:len(y_range))
ϵ_y_dict = Dict(ϵ_y_range[i] => i for i = 1:len(ϵ_y_range))
range_dicts = (π = π_dict, y=y_dict,ϵ_y_lead=ϵ_y_dict)
ranges = (π_lag =π_range, y_lag = y_range, ϵ_y_range = ϵ_y_range)

n_states=1
for vv in 1:len(ranges)
	global n_states *=len(ranges[vv])
end
n_states = len(y_range)*len(ϵ_y_range)

state_array = Array(reshape(Array(1:n_states),(len(ϵ_y_range),len(y_range))))

state_t_vars = [:y_lag, :ϵ_y]
state_tp1_vars = [:y, :ϵ_y_lead]
endog = indices.endognames
endog_lead = Symbol.(String.(indices.endognames) .* "_lead")
endog_lag = Symbol.(String.(indices.endognames) .* "_lag")

for it in 1:100
    display(join(["Running through iteration ", it]))
	global df_grid = create_df_grid(state_grid, beliefs, indices)
	n_states = nrow(df_grid)
	all_states = 1:n_states
	df_grid.state_t = all_states
	df_grid.state_tp1 = Int.(zeros(nrow(df_grid)))
	df_grid.uncon_prob = zeros(nrow(df_grid))
	@time df_grid = step_1_map(df_grid, beliefs, indices, options, sval = 0.0)

	# Make sure that the step1 variables correspond to grid points
	for vv in 1:len(endog)
		df_grid[:,endog[vv]] = map(x -> findnearest(ranges[endog_lag[vv]],x), Array(df_grid[:,endog[vv]]))
	end
	# Find step 2 points and transition probabilities
	@time global df_grid_new, loss_weights, probs_df = create_trans_probs(df_grid, transition_probabilities,
		indices, shock_range, parallel = true)

	# Classify state for t+1
	# Number of rows in df_grid is number of possible states
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
	global CP = create_kernel(df_grid_new, n_states)

	#display(heatmap(CP))
	CP ./= sum(CP,dims = 2)
	MP = MC_stationary_fast(CP)
	#@elapsed MP = MC_stationary(CP)
	MP = MP./sum(MP)
	plot(1:length(MP),MP, title = "Marginal prob of states")

	# Use Bayes Rule to find unconditional probabilities
	MP_mat = repeat(MP, outer = [1, len(MP)])
	UP = MP_mat.*CP
	UP ./= sum(UP)

	# Extract unconditional probability associated with each (t,t+1) pair
	for st in 1:nrow(df_grid_new)
		df_grid_new[st,:uncon_prob] = UP[df_grid_new.state_t[st],df_grid_new.state_tp1[st]]
	end

	display(plot(df_grid_new.y, df_grid_new.uncon_prob, title = "Unconditional prob of y"))
	df_train = df_grid_new[(df_grid_new.uncon_prob .> 1e-16),:]

	# Identify the remaining states
	states = unique(df_train.state_t)

	loss_weights = Matrix(transpose(df_train.uncon_prob))
	loss_weights ./=(sum(loss_weights)/size(loss_weights,2))
	loss_weights = repeat(loss_weights, len(indices.expectnames_all))

	@time global beliefs = update_beliefs(df_train, beliefs, indices, options,
		epochs = 10, cycles = 5, verbose = true, weights = loss_weights)
	if true
		starting_point = deepcopy(upper)
		starting_point[:ϵ_y] = 0.0
		paths1 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	    	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	    	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
		starting_point = deepcopy(lower)
		starting_point[:ϵ_y] = 0.0
		paths2 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	    	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	    	plot_vars = (:p, :y, :Ep_lead), y_lim = [-3.0,3.0])
		plot(layout = (2,1), title = join(["Iteration: ", it]))
		plot!(paths1.π, subplot =1)
		plot!(paths1.y, subplot =1)
		plot!(paths2.π, subplot =2)
		display(plot!(paths2.y, subplot =2))
	end
end
#save("networks/meandyn_2eq_24nodes.jld2", "beliefs", beliefs)
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
