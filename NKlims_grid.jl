cd("/Users/julianashwin/Documents/GitHub/EcoNNet.jl")

using Distributed
#addprocs(4)
#rmprocs(4)
workers()
@everywhere include("src/EcoNNet_dev.jl")

"""
Define some options.
Expectations are denoted by a capital E prefix,
leads and lags are specified as a `_lead` or `_lag` suffix
"""

@everywhere options = EcoNNetOptions(infoset = [:y_lag, :ϵ_y],
    expectations = [:Eπ,:Ey,:Eπ_lead, :Ey_lead],
    endogenous = [:π, :y],
    exogenous = [:ϵ_y],
    states = [:y_lag, :ϵ_y],
	auxiliary = [:r],
    N = 500000, num_nodes = 24, activation = σ, window = 40000)

@everywhere beliefs = initialise_beliefs(options)

"""
Define the parameters as a Named Tuple.
"""

@everywhere par = (β = 0.95, κ = 0.1,
	η = 0.5, σ = 0.25, ϕ_π = 1.5, ϕ_y = 0.5,
	ρ_π = 0.5, σ_π = 0.001, ρ_y = 0.5, σ_y = 0.005,
	R_lim = -0.02, π_lim = -0.03, y_lim = -0.05);

"""
State the equilibrium conditions of the model as a function which returns
    a vector of zeros
"""

@everywhere function equilibrium_conditions_fast(F::Array{Float64,1},
    x::Array{Float64,1},states_input::Array{Float64,1},predictions_input::Array{Float64,1})
    # Manually unpack the states
    #π_lag::Float64 = states_input[1]
    y_lag::Float64 = states_input[1]
    #ϵ_π::Float64 = states_input[3]
	ϵ_y::Float64 = states_input[2]
    # and the predictions
    #Ep::Float64 = predictions_input[1]
    #Ey::Float64 = predictions_input[2]
    Eπ_lead::Float64 = predictions_input[3]
	Ey_lead::Float64 = predictions_input[4]
    # and the endogenous variables
    π::Float64 = x[1]
    y::Float64 = x[2]

    # p[t] = ϕ_pp*pe[t+1] + ϕ_py*y[t] + α_2*y[t]^2 - α_3*y[t]^3 + ϵ_p[t]
    F[1] =  max(par.β*Eπ_lead + par.κ*y,par.π_lim) - π ::Float64
    # y[t] = η*y[t-1] + ϕ_yp*p[t] + ϵ_y[t]
	r = max(par.ϕ_π*π + par.ϕ_y*y, par.R_lim)
    F[2] = max((1-par.η)*Ey_lead + par.η*y_lag - par.σ*(r - Eπ_lead) , par.y_lim)+ ϵ_y - y ::Float64

    return F
end


"""
Define equations of motion under perfect foresight as a useful comparison
	Input - first lagged state variables, then endogenous variables
"""
function perfect_foresight(inputs)
    # Manually unpack the states (endogenous variables in same order as in options)
	# Lagged state variables first, then current condogenous
    π_lag = max(inputs[1], par.π_lim)
    y_lag = max(inputs[2], par.y_lim)
	π_t = max(inputs[3], par.π_lim)
	y_t = max(inputs[4], par.y_lim)



    # π[t+1] = 1/β*(π[t] - κ*y[t])
    π_new = max((1/par.β)*(π_t - par.κ*y_t),par.π_lim)
    # y[t+1] = η*y[t] - σ*(ϕ_π*π[t+1] + α π[t+1]^3 - π[t+2])
	r_t = max(par.ϕ_π*π_t + par.ϕ_y*y_t, par.R_lim)
	y_new = (1/(1-par.η))*(y_t - par.η*y_lag + par.σ*(
		r_t - π_new))
	y_new = max(y_new,par.y_lim)
	# Impose upper and lower bounds to allow plotting
	π_new = min(π_new,max(π_new,-1e6),1e6)
	y_new = min(y_new,max(y_new,-1e6),1e6)

    outputs = [π_new,y_new]

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
function steady_states(F::Array{Float64,1},x::Array{Float64,1})
    π::Float64 = x[1]
    y::Float64 = x[2]
    F[1]::Float64 = max(par.β*π + par.κ*y,par.π_lim) - π
    F[2]::Float64 = max((1-par.η)*y + par.η*y - par.σ*(
		max(par.ϕ_π*π + par.ϕ_y*y, par.R_lim) - π),par.y_lim) - y
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
sstates = nlsolve(steady_states, [-0.02, -0.02])
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
function ZLB_bind_condition(y)
	π = (par.R_lim - par.ϕ_y*y)/par.ϕ_π
end
function NKPC_condition(y)
    π = max((par.κ/(1-par.β))*y, par.π_lim)
end
function EE_normal_condition(y)
	π = (par.ϕ_y/(1-par.ϕ_π))*y
end
function EE_ZLB_condition(y)
	π = par.R_lim
end
function EE_ZLB_condition(y)
	π = par.R_lim
end
EE_kink = par.R_lim*(1-par.ϕ_π)/(par.ϕ_π*par.ϕ_y+(1-par.ϕ_π)*par.ϕ_y)

## Steady state conditions
gr()
plot_points = -1:0.001:1
ss_plot = plot(xlabel = "Output", xlims = (-0.07,0.07),
    ylabel = "Inflation", ylims = (-0.07, 0.07),legend=:topright)
plot!(plot_points, ZLB_bind_condition.(plot_points),
	fillrange=[minimum(plot_points)*ones(len(plot_points))], fillalpha = 0.5,
	 color = :paleturquoise1, label = "ZLB")
plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
plot!(par.y_lim:0.001:EE_kink,EE_normal_condition.(par.y_lim:0.001:EE_kink),
 	label = "Euler Equation", color = :green)
plot!(par.y_lim:0.001:EE_kink,EE_ZLB_condition.(par.y_lim:0.001:EE_kink),
 	label = "", color = :green)
plot!([par.y_lim,par.y_lim],[par.R_lim,-0.07],
 	label = "", color = :green)
display(ss_plot)

## Plot perfect foresight paths
phase_plot = ss_plot
initial_ss = deepcopy(central)
starts = [(π=0.003,y=0.003,periods=5,arrows=[4]),
	(π=-0.003,y=-0.003,periods=5,arrows=[4]),
	(π=0.01,y=-0.003,periods=7,arrows=[6]),
	(π=-0.01,y=0.003,periods=7,arrows=[6]),
	#(π=-0.025,y=-0.03,periods=3,arrows=[2]),
	(π=-0.045,y=-0.03,periods=6,arrows=[5]),
	(π=-0.03,y=-0.05,periods=12,arrows=[7,11]),
	(π=-0.0292,y=-0.04,periods=7,arrows=[6]),
	(π=-0.015,y=0.005,periods=6,arrows=[5]),
	(π=-0.0175,y=0.005,periods=6,arrows=[5])]
for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods])
	phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
		v_points = 1:start[:periods])
end
display(phase_plot)
#savefig("figures/phase_NKlim_perf.png")

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
# Initialise beliefs
@everywhere beliefs = initialise_beliefs(options)
s = initialise_df(s, lower, gap = 500, steadystate_alt = upper)
s.ϵ_y = simulate_ar(par.ρ_y, par.σ_y, options.N, Normal())
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss, cutoff = false)

starting_point = deepcopy(upper)
starting_point[:ϵ_y] = 0.0
paths1 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	plot_vars = (:p, :y, :Ep_lead), y_lim = [-0.08,0.04])
starting_point = deepcopy(lower)
starting_point[:ϵ_y] = 0.0
paths2 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	plot_vars = (:p, :y, :Ep_lead), y_lim = [-0.08,0.04])
plot(layout = (2,1))
plot!(paths1.π, subplot =1)
plot!(paths1.y, subplot =1)
plot!(paths2.π, subplot =2)
display(plot!(paths2.y, subplot =2))

"""
Use Network in grid solution method
"""
@everywhere par = (β = 0.95, κ = 0.1,
	η = 0.5, σ = 0.25, ϕ_π = 1.5, ϕ_y = 0.5,
	ρ_π = 0.5, σ_π = 0.001, ρ_y = 0.5, σ_y = 0.01,
	R_lim = -0.02, π_lim = -0.03, y_lim = -0.05);


p_y = (1 + par.ρ_y)/2
q_y = p_y
N = 21
ψ_y = sqrt(par.σ_y^2/(1 - par.ρ_y^2))*sqrt(N-1)
ϵ_y_range, ϵ_y_kernel = Rouwenhorst(p_y,q_y,ψ_y,N)
gr()
heatmap(ϵ_y_range,
    ϵ_y_range, ϵ_y_kernel,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="shock(t+1)", ylabel="shock(t)",
    title="My title")


display(join(["Order in grid should be ",indices.statenames_all]))
π_range = Array(-0.04:0.001:0.04)
y_range = Array(-0.08:0.001:0.04)
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
Grid-based Mean Dynamics
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
		endog_def = df_grid[st,endog]
		sel_rows = BitArray(((df_grid_new.y .== state_def.y_lag).*
			(df_grid_new.ϵ_y_lead .== state_def.ϵ_y)))
		df_grid_new[sel_rows,:state_tp1] .= st
		df_grid_new[sel_rows,:π_lead] .= endog_def.π
		df_grid_new[sel_rows,:y_lead] .= endog_def.y
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
	df_grid_new.uncon_prob = len(ϵ_y_range).*df_grid_new.uncon_prob
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
	    	plot_vars = (:p, :y, :Ep_lead), y_lim = [-0.08,0.04])
		starting_point = deepcopy(lower)
		starting_point[:ϵ_y] = 0.0
		paths2 = irf(:ϵ_y, starting_point, beliefs, shock_period = 10, periods = 100,
	    	magnitude = 0.0, persistence = par.ρ_y, show_plot = false,
	    	plot_vars = (:p, :y, :Ep_lead), y_lim = [-0.08,0.04])
		plot(layout = (2,1), title = join(["Iteration: ", it]))
		plot!(paths1.π, subplot =1)
		plot!(paths1.y, subplot =1)
		plot!(paths2.π, subplot =2)
		display(plot!(paths2.y, subplot =2))
	end
end















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

paths1 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 100,
	magnitude = -0.0, persistence = par.ρ_y, show_plot = true,plot_vars=[:y,:π])
paths2 = irf(:ϵ_y, lower_stoch, beliefs, shock_period = 5, periods = 100,
	magnitude = -0.01, persistence = par.ρ_y, show_plot = true,
	plot_vars=[:y,:π,:Ey, :Eπ])

"""
Plot IRFs
"""

paths1 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 100,
	magnitude = -0.005, persistence = par.ρ_y, show_plot = false)
paths2 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 100,
	magnitude = 0.005, persistence = par.ρ_y, show_plot = false)
paths3 = irf(:ϵ_y, upper_stoch, beliefs, shock_period = 5, periods = 100,
	magnitude = -0.03, persistence = par.ρ_y, show_plot = false)
paths4 = irf(:ϵ_y, lower_stoch, beliefs, shock_period = 5, periods = 100,
	magnitude = 0.03, persistence = par.ρ_y, show_plot = false)
plot(paths1.π, label ="Inflation", xlabel = "Periods", legend = :bottomright,ylims = (-0.08,0.08))
plot!(paths1.y, label ="Output");plot!(paths1.ϵ_y, label ="Shock")
plot(paths2.π, label ="Inflation", xlabel = "Periods", legend = false,ylims = (-0.08,0.08))
plot!(paths2.y, label ="Output");plot!(paths2.ϵ_y, label ="Shock")
plot(paths3.π, label ="Inflation", xlabel = "Periods", legend = false,ylims = (-0.08,0.08))
plot!(paths3.y, label ="Output");plot!(paths3.ϵ_y, label ="Shock")
plot(paths4.π, label ="Inflation", xlabel = "Periods", legend = false,ylims = (-0.08,0.08))
plot!(paths4.y, label ="Output");plot!(paths4.ϵ_y, label ="Shock")

# All plots together
irf_plot = plot(layout = (2,2),ylims = (-0.06,0.06),size=(1000,600))
plot!(paths1.π, label ="Inflation", xlabel = "Periods", legend = :topright, subplot=1)
plot!(paths1.y, label ="Output");plot!(paths1.ϵ_y, label ="Shock")
plot!(paths2.π, label ="Inflation", xlabel = "Periods", legend = false, subplot=2)
plot!(paths2.y, label ="Output", subplot=2);plot!(paths2.ϵ_y, label ="Shock", subplot=2)
plot!(paths3.π, label ="Inflation", xlabel = "Periods", legend = false, subplot=3)
plot!(paths3.y, label ="Output", subplot=3);plot!(paths3.ϵ_y, label ="Shock", subplot=3)
plot!(paths4.π, label ="Inflation", xlabel = "Periods", legend = false, subplot=4)
plot!(paths4.y, label ="Output", subplot=4);plot!(paths4.ϵ_y, label ="Shock", subplot=4)

display(irf_plot)
savefig("figures/NK_irf_sim.png")


"""
Plot phase diagram
"""

gr()
plot_points = -1:0.001:1
ss_plot = plot(xlabel = "Output", xlims = (-0.07,0.07),
    ylabel = "Inflation", ylims = (-0.07, 0.07),legend=:topright)
plot!(plot_points, ZLB_bind_condition.(plot_points),
	fillrange=[minimum(plot_points)*ones(len(plot_points))], fillalpha = 0.5,
	 color = :paleturquoise1, label = "ZLB")
plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
plot!(par.y_lim:0.001:EE_kink,EE_normal_condition.(par.y_lim:0.001:EE_kink),
 	label = "Euler Equation", color = :green)
plot!(par.y_lim:0.001:EE_kink,EE_ZLB_condition.(par.y_lim:0.001:EE_kink),
 	label = "", color = :green)
plot!([par.y_lim,par.y_lim],[par.R_lim,-0.07],
 	label = "", color = :green)
display(ss_plot)

## Plot perfect foresight paths
phase_plot = ss_plot
initial_ss = deepcopy(central)
starts = [(π=-0.04,y=-0.06,ϵ_y=-0.0,periods=2,arrows=[1]),
	(π=-0.04,y=-0.018,ϵ_y=0.0,periods=43,arrows=[27,32]),
	(π=-0.04,y=-0.016,ϵ_y=0.0,periods=43,arrows=[]),
	(π=0.035,y=0.02,ϵ_y=0.0,periods=12,arrows=[6]),
	(π=0.0,y=-0.039,ϵ_y=0.0,periods=60,arrows=[58]),
	(π=0.0,y=-0.037,ϵ_y=0.0,periods=100,arrows=[18,90]),
	(π=-0.0175,y=0.05,ϵ_y=0.0,periods=60,arrows=[35])
	]
for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = 100,
		magnitude = start[:ϵ_y], persistence = par.ρ_y, show_plot = false)
	phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
		v_points = 1:start[:periods])
end
display(phase_plot)
savefig("figures/NK_phase_sim.png")


"""
Compare PLM and ALM
"""
s.Eπ_error = s.Eπ-s.π
s.Ey_error = s.Ey-s.y
s.Eπ_lead_error = vcat((s.Eπ_lead[1:(options.N-1)]-s.π[2:(options.N)]),[0.])
plot(s.Eπ_lead_error)

plot_range=75001:100000
y_range = Array(-0.1:0.0002:0.04)
π_range= Array(-0.04:0.0002:0.02)
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
savefig("figures/NK_heatmap_sim.png")


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
