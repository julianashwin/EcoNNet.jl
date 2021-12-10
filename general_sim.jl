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

@everywhere options = EcoNNetOptions(infoset = [:p_lag, :y_lag, :ϵ_p, :ϵ_y],
    expectations = [:Ep,:Ey,:Ep_lead],
    endogenous = [:p, :y],
    exogenous = [:ϵ_p,:ϵ_y],
    states = [:y_lag, :ϵ_p, :ϵ_y],
	auxiliary = [],
	burnin = 100000,
    N = 1000000, num_nodes = 24, activation = relu, max_iter = 20,
	burnin_use_net = false, window = 99999, learning_gap = 100000,
	plotting_gap = 100000, plot_vars = [:p, :y, :Ep, :Ey])

@everywhere beliefs = initialise_beliefs(options)

"""
Define the parameters as a Named Tuple.
"""

@everywhere par = (ϕ_pp = 0.95, ϕ_py = 0.5,
	ϕ_yy = 0.9, ϕ_yp = 0.1,
    α_2 = 0.0, α_3 = 0.1125,
	ρ_y = 0.5, σ_y = 0.3, ρ_p = 0.5, σ_p = 0.3);


"""
State the equilibrium conditions of the model as a function which returns
    a vector of zeros
"""

@everywhere function equilibrium_conditions_fast(F::Array{Float64,1},
    x::Array{Float64,1},states_input::Array{Float64,1},predictions_input::Array{Float64,1})
    # Manually unpack the states
    #p_lag::Float64 = states_input[1]
    y_lag::Float64 = states_input[1]
    ϵ_p::Float64 = states_input[2]
	ϵ_y::Float64 = states_input[3]
    # and the predictions
    #Ep::Float64 = predictions_input[1]
    #Ey::Float64 = predictions_input[2]
    Ep_lead::Float64 = predictions_input[3]
    # and the endogenous variables
    p::Float64 = x[1]
    y::Float64 = x[2]

    # p[t] = ϕ_pp*pe[t+1] + ϕ_py*y[t] + α_2*y[t]^2 - α_3*y[t]^3 + ϵ_p[t]
    F[1] =  par.ϕ_pp*Ep_lead + par.ϕ_py*y - par.α_3*y^3 + ϵ_p - p ::Float64
    # y[t] = η*y[t-1] + ϕ_yp*p[t] + ϵ_y[t]
    F[2] = par.ϕ_yy*y_lag + par.ϕ_yp*p + ϵ_y - y ::Float64

    return F
end


"""
Define equations of motion under perfect foresight as a useful comparison
"""
function perfect_foresight(inputs)
    # Manually unpack the states (endogenous variables in same order as in options)
    p_lag = inputs[1]
    y_lag = inputs[2]
    # and the predictions
    # π[t+1] = 1/β*(π[t] - κ*y[t])
    p = (1/par.ϕ_pp)*(p_lag - par.ϕ_py*y_lag + par.α_3*y_lag^3)
    # y[t+1] = η*y[t] - σ*(ϕ_π*π[t+1] + α π[t+1]^3 - π[t+2])
	y = par.ϕ_yy*y_lag + par.ϕ_yp*p
	# Impose upper and lower bounds to allow plotting
	p = min(p,max(p,-1e6),1e6)
	y = min(y,max(y,-1e6),1e6)

    outputs = [p,y]

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
function p_condition(y)
    p = 1/(1-par.ϕ_pp)*(par.ϕ_py*y - par.α_3*y^3)
end
function y_condition(p)
    y = (par.ϕ_yp/(1 - par.ϕ_yy))*p
end
function steady_states(F::Array{Float64,1},x::Array{Float64,1})
    p::Float64 = x[1]
    y::Float64 = x[2]
    F[1]::Float64 = par.ϕ_pp*p + par.ϕ_py*y - par.α_3*y^3 - p
    F[2]::Float64 = par.ϕ_yy*y + par.ϕ_yp*p - y
    return F
end

# Exogenous processes
ss = Dict{Symbol,Float64}()
ss[:ϵ_p] = 0.0
ss[:ϵ_y] = 0.0
# upper steady state
sstates = nlsolve(steady_states, [2.0, 2.0])
upper = Dict{Symbol,Float64}();
upper[:p] = sstates.zero[1];
upper[:y] = sstates.zero[2];
upper[:Ep_lead] = upper[:p];
# central steady state
sstates = nlsolve(steady_states, [0.0, 0.0])
central = Dict{Symbol,Float64}();
central[:p] = sstates.zero[1];
central[:y] = sstates.zero[2];
central[:Ep_lead] = central[:p];
# lower steady state
lower = Dict{Symbol,Float64}();
sstates = nlsolve(steady_states, [-2.0, -2.0])
lower[:p] = sstates.zero[1];
lower[:y] = sstates.zero[2];
lower[:Ep_lead] = lower[:p];

s = initialise_df(s, upper);
s = initialise_df(s, ss);
s = initialise_df(s, lower, gap = 2, steadystate_alt = upper)


"""
Plot steady state conditions and perfect foresight paths
"""

#using Plotly, PlotlyJS
pyplot()
plot_points = -4.0:0.01:4.0;
function plot_ss(plot_points)
	ss_plot = plot(xlabel = L"y_{t-1}", xlims = (minimum(plot_points),maximum(plot_points)),
    	ylabel = L"p_t", ylims = (-8.5,8.5),legend=:bottomright, yguidefontrotation=-90)
		plot!(plot_points,p_condition.(plot_points), label = "p condition", color = :black)
		display(plot!(y_condition.(plot_points),plot_points, label = "y condition", color = :green))
		# Plot perfect foresight paths
	return ss_plot
end

starts = [(p=-1.5,y=3.0,periods=100,arrows=[4,60]),
	(p=-1.808,y=3.0,periods=12,arrows=[2]),
	(p=-2.5,y=3.0,periods=62,arrows=[10,60]),
	(p=-3.75,y=3.0,periods=62,arrows=[10,60]),
	(p=-3.338,y=3.0,periods=25,arrows=[6,12]),
	(p=1.5,y=-3.0,periods=100,arrows=[4,60]),
	(p=1.808,y=-3.0,periods=12,arrows=[2]),
	(p=2.5,y=-3.0,periods=62,arrows=[10,60]),
	(p=3.75,y=-3.0,periods=62,arrows=[10,60]),
	(p=3.338,y=-3.0,periods=25,arrows=[10,12])]

#@everywhere par = (ϕ_pp = 1.3968, ϕ_py = -0.3069,
#	ϕ_yy = 0.9787, ϕ_yp = 0.0201,
#    α_2 = 0.0, α_3 = 0.0,
#	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);
#@everywhere par = (ϕ_pp = 1.0212, ϕ_py = -0.0075,
#	ϕ_yy = 1.0183, ϕ_yp = -0.0423,
#    α_2 = 0.0, α_3 = 0.0,
#	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);
#starts = [(p=5.0,y=3.0,periods=100,arrows=[10]),
#	(p=-5.0,y=3.0,periods=300,arrows=[10]),
#	(p=5.0,y=-3.0,periods=100,arrows=[10]),
#	(p=-5.0,y=-3.0,periods=300,arrows=[10])]


#for tt in 1:50
plot_ss(plot_points)
initial_ss = deepcopy(central)
for start in starts
	initial_ss[:p] = start[:p]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods])

	if start == starts[1]
		phase_arrow_plot(paths1, [:y,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 1:(start[:periods]-1), label = "Perfect foresight paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths1, [:y,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 1:(start[:periods]-1), label = "", arrow_size = .5,
			final_arrow = true)
	end
end
plot!(size = (600,400))
#plot!(title = "Perfect foresight")
savefig("figures/general/phase_compsink_perf.pdf")



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

noise_p = par.σ_p*randn(nrow(s))
s.ϵ_p = simulate_ar(par.ρ_p, par.σ_p, options.N, noise_p)
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
options.plot_vars = [:p, :y, :Ep, :Ey]

# Simulate the learning for a set number of periods
noise_p = par.σ_p*randn((options.N - options.burnin))
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_p[(options.burnin+1):options.N] = simulate_ar(par.ρ_p, par.σ_p, options.N - options.burnin, noise_p)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs, s = simulate_learning((options.burnin+1):options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-5999):(500000-4999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.p[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))
#savefig("figures/general/general_sim_series.pdf")

export_df = s[options.N-99999:options.N,:]
rename!(export_df, replace.(names(export_df), "ϵ" => "epsilon"))

#CSV.write("estimation/general/general_sim.csv", export_df)


"""
Plot phase diagram
"""

pyplot()
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(p=-1.5,y=3.0,periods=100,arrows=[4,60]),
	(p=1.5,y=-3.0,periods=100,arrows=[4,60]),
	(p=0.1,y=0.1,periods=100,arrows=[25]),
	(p=-0.2,y=-0.2,periods=100,arrows=[25]),
	]
for start in starts
	initial_ss[:p] = start[:p]; initial_ss[:y] = start[:y];
	paths = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = start[:periods],
		magnitude = 0.0, persistence = par.ρ_y, show_plot = false)
	paths.y_lag = cat(start[:y],paths.y[1:(start.periods -1 )],dims=1)
	if start == starts[1]
		phase_arrow_plot(paths, [:y_lag,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 2:(start[:periods]), label = "Equilibrium paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths, [:y_lag,:p], arrow_points=start[:arrows], h_points = 2:(start[:periods]),
			v_points = 2:start[:periods], arrow_size = .5, final_arrow = true)
		end
end
plot!(size = (600,400))
savefig("figures/general/phase_compsink_sim.pdf")



















"""
Real sink
"""

@everywhere par = (ϕ_pp = 0.95, ϕ_py = 0.5,
	ϕ_yy = 0.6, ϕ_yp = 0.1,
    α_2 = 0.0, α_3 = 0.1125,
	ρ_y = 0.0, σ_y = 0.2, ρ_p = 0.95, σ_p = 0.05);

ystar = 0.0
a_11 = (1 - (par.ϕ_py - 3*ystar^2*par.α_3)*par.ϕ_yp)/par.ϕ_pp
a_12 = (-(par.ϕ_py - 3*ystar^2*par.α_3)*par.ϕ_yy)/par.ϕ_pp
a_21 = par.ϕ_yp
a_22 = par.ϕ_yy
A_matrix = [a_11 a_12;
	a_21 a_22]
eigen(A_matrix)

sstates = nlsolve(steady_states, [2.0, 2.0])
upper[:p] = sstates.zero[1];
upper[:y] = sstates.zero[2];
upper[:Ep_lead] = upper[:p];
# lower steady state
sstates = nlsolve(steady_states, [-2.0, -2.0])
lower = Dict{Symbol,Float64}();
lower[:p] = sstates.zero[1];
lower[:y] = sstates.zero[2];
lower[:Ep_lead] = lower[:p];


"""
Perfect foresight phase
"""

pyplot()
plot_points = -4.0:0.01:4.0;
function plot_ss(plot_points)
	ss_plot = plot(xlabel = L"y_{t-1}", xlims = (minimum(plot_points),maximum(plot_points)),
    	ylabel = L"p_t", ylims = (-10.0,10.0),legend=:bottomright, yguidefontrotation=-90)
		plot!(plot_points,p_condition.(plot_points), label = "p condition", color = :black)
		display(plot!(y_condition.(3.5.*plot_points),3.5.*plot_points, label = "y condition", color = :green))
		# Plot perfect foresight paths
	return ss_plot
end

starts = [(p=3.0,y=1.0,periods=20,arrows=[1]),
	(p=-3.0,y=-1.0,periods=20,arrows=[1]),
	(p=2.0,y=1.0,periods=20,arrows=[1,4]),
	(p=-2.0,y=-1.0,periods=20,arrows=[1,4]),
	(p=-2.5,y=3.75,periods=18,arrows=[1]),
	(p=2.5,y=-3.75,periods=18,arrows=[1]),
	(p=-2.5,y=3.85,periods=30,arrows=[1]),
	(p=2.5,y=-3.85,periods=30,arrows=[1]),
	#(p=7.0,y=-3.5,periods=40,arrows=[4]),
	#(p=-7.0,y=3.5,periods=40,arrows=[4]),
	#(p=-8.32,y=3.5,periods=10,arrows=[3]),
	#(p=0.0,y=3.0,periods=10,arrows=[3])
	]
starts = [(p=-1.0,y=4.1,periods=20,arrows=[4]),
	(p=1.0,y=-4.1,periods=20,arrows=[4]),
	(p=-1.0,y=3.8,periods=18,arrows=[4]),
	(p=1.0,y=-3.8,periods=18,arrows=[4]),
	(p=6.5,y=-2.0,periods=30,arrows=[4]),
	(p=-6.5,y=2.0,periods=30,arrows=[4]),
	(p=7.0,y=-3.5,periods=40,arrows=[4]),
	(p=-7.0,y=3.5,periods=40,arrows=[4]),
	(p=1.55,y=3.5,periods=6,arrows=[2]),
	(p=-1.55,y=-3.5,periods=6,arrows=[2]),
	(p=8.32,y=-3.5,periods=10,arrows=[3]),
	(p=-8.32,y=3.5,periods=10,arrows=[3]),
	]
plot_ss(plot_points)
initial_ss = deepcopy(central)
for start in starts
	initial_ss[:p] = start[:p]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods])

	if start == starts[1]
		phase_arrow_plot(paths1, [:y,:p], arrow_points=start[:arrows], h_points = 1:start[:periods],
			v_points = 1:(start[:periods]), label = "Perfect foresight paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths1, [:y,:p], arrow_points=start[:arrows], h_points = 1:start[:periods],
			v_points = 1:(start[:periods]), label = "", arrow_size = .5,
			final_arrow = true)
	end
end
plot!(size = (600,400))
savefig("figures/general/phase_realsink_perf.pdf")


"""
Re-train beliefs
"""
# Initialise beliefs
@everywhere beliefs = initialise_beliefs(options)
s = initialise_df(s, central, gap = 500, steadystate_alt = central)
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss)

# Simulate the learning for a set number of periods
noise_p = par.σ_p*randn(nrow(s)) + alternating_shocks(nrow(s), gap=400, mag = 1.0)
s.ϵ_p = simulate_ar(par.ρ_p, par.σ_p, options.N, noise_p)
plot(s.ϵ_p[1:1000])
noise_y = par.σ_y*randn(nrow(s))
s.ϵ_y = simulate_ar(par.ρ_y, par.σ_y, options.N, noise_y)
plot(s.ϵ_y[1:1000])
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], upper, gap = 500, steadystate_alt = upper);
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
@time beliefs, s = simulate_learning((options.burnin+1):options.N, s, beliefs, indices, options)

# Repeat learning with new shocks
noise_p = par.σ_p*randn((options.N - options.burnin)) +
	alternating_shocks((options.N - options.burnin), gap=300, mag = 1.0)
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_p[(options.burnin+1):options.N] = simulate_ar(par.ρ_p, par.σ_p, options.N - options.burnin, noise_p)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs, s = simulate_learning((options.burnin+1):options.N, s, beliefs, indices, options)


# Plot simulated time series
pyplot()
plot_range = (500000-5999):(500000-4999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.p[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))

# Check IRFs
paths1 = irf(:ϵ_y, upper, beliefs, shock_period = 5, periods = 100,
	magnitude = -0.5, persistence = par.ρ_y, show_plot = false)
paths2 = irf(:ϵ_y, lower, beliefs, shock_period = 5, periods = 100,
	magnitude = 0.5, persistence = par.ρ_y, show_plot = false)
plot(paths1.p, label =L"p_t", xlabel = "Periods", legend = :bottomright)
plot!(paths1.y, label =L"y_t");
plot!(paths1.ϵ_y, label =L"\epsilon_{y,t}")
plot(paths2.p, label =L"p_t", xlabel = "Periods", legend = false)
plot!(paths2.y, label =L"y_t");
plot!(paths2.ϵ_y, label =L"\epsilon_{y,t}")


"""
Plot phase diagram
"""
# Solution phase diagram
pyplot()
plot_ss(plot_points)
#plot!(ylims = (-12.0,12.0))
initial_ss = deepcopy(central)
starts = [#(p=0.0,y=2.0,periods=15,arrows=[2]),
	(p=5.0,y=3.5,periods=15,arrows=[2]),
	(p=-5.0,y=-3.5,periods=15,arrows=[2]),
	#(p=-3.0,y=0.0,periods=15,arrows=[2]),
	#(p=0.0,y=-1.0,periods=100,arrows=[25]),
	#(p=5.0,y=0.0,periods=100,arrows=[25]),
	#(p=-5.0,y=-0.0,periods=100,arrows=[25]),
	]
for start in starts
	initial_ss[:p] = start[:p]; initial_ss[:y] = start[:y];
	paths = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = start[:periods],
		magnitude = 1.0, persistence = par.ρ_y, show_plot = false)
	paths.y_lag = cat(start[:y],paths.y[1:(start.periods -1 )],dims=1)
	if start == starts[1]
		phase_arrow_plot(paths, [:y_lag,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 2:(start[:periods]), label = "Equilibrium paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths, [:y_lag,:p], arrow_points=start[:arrows], h_points = 2:(start[:periods]),
			v_points = 2:start[:periods], arrow_size = .5, final_arrow = true)
		end
end
plot!(size = (600,400))
savefig("figures/general/phase_realsink_sim.pdf")














"""
Complex source
"""

@everywhere par = (ϕ_pp = 0.95, ϕ_py = 0.5,
	ϕ_yy = 1.5, ϕ_yp = 0.1,
    α_2 = 0.0, α_3 = 0.1125,
	ρ_y = 0.5, σ_y = 0.5, ρ_p = 0.5, σ_p = 1.5);

a_11 = (1 - par.ϕ_py*par.ϕ_yp)/par.ϕ_pp
a_12 = - (par.ϕ_py*par.ϕ_yy)/par.ϕ_pp
a_21 = par.ϕ_yp
a_22 = par.ϕ_yy
A_matrix = [a_11 a_12;
	a_21 a_22]
eigen(A_matrix)

sstates = nlsolve(steady_states, [2.0, 2.0])
upper[:p] = sstates.zero[1];
upper[:y] = sstates.zero[2];
upper[:Ep_lead] = upper[:p];
# lower steady state
lower = Dict{Symbol,Float64}();
lower[:p] = sstates.zero[1];
lower[:y] = sstates.zero[2];
lower[:Ep_lead] = lower[:p];


"""
Perfect foresight phase
"""
pyplot()
plot_points = -4.0:0.01:4.0;
function plot_ss(plot_points)
	ss_plot = plot(xlabel = L"y_{t-1}", xlims = (minimum(plot_points),maximum(plot_points)),
    	ylabel = L"p_t", ylims = (-16.0,16.0),legend=:bottomright, yguidefontrotation=-90)
		plot!(plot_points,p_condition.(plot_points), label = "p condition", color = :black)
		display(plot!(y_condition.(5.0.*plot_points),5.0.*plot_points, label = "y condition", color = :green))
		# Plot perfect foresight paths
	return ss_plot
end

starts = [(p=0.05,y=0.05,periods=20,arrows=[6]),
	(p=-0.05,y=-0.05,periods=20,arrows=[6]),
	(p=0.1,y=-0.05,periods=20,arrows=[16]),
	(p=-0.1,y=0.05,periods=20,arrows=[16]),
	(p=1.5,y=-0.5755804,periods=25,arrows=[4]),
	(p=-1.5,y=0.5755804,periods=25,arrows=[4]),
	(p=1.5,y=-0.577,periods=17,arrows=[4]),
	(p=-1.5,y=0.577,periods=17,arrows=[4]),
	(p=1.5,y=-0.574,periods=18,arrows=[4]),
	(p=-1.5,y=0.574,periods=18,arrows=[4]),
	]
plot_ss(plot_points)
plot!(legend = :topright)
initial_ss = deepcopy(central)
for start in starts
	initial_ss[:p] = start[:p]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods])

	if start == starts[1]
		phase_arrow_plot(paths1, [:y,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 1:(start[:periods]-1), label = "Perfect foresight paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths1, [:y,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 1:(start[:periods]-1), label = "", arrow_size = .5,
			final_arrow = true)
	end
end
plot!(size = (600,400))
savefig("figures/general/phase_compsource_perf.pdf")


"""
Re-train beliefs
"""
# Initialise beliefs
@everywhere beliefs = initialise_beliefs(options)
s = initialise_df(s, lower, gap = 500, steadystate_alt = upper)
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss)

# Simulate the learning for a set number of periods
noise_p = par.σ_p*randn(nrow(s))
s.ϵ_p = simulate_ar(par.ρ_p, par.σ_p, options.N, noise_p)
noise_y = par.σ_y*randn(nrow(s))
s.ϵ_y = simulate_ar(par.ρ_y, par.σ_y, options.N, noise_y)
plot(s.ϵ_y[1:1000])
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], lower, gap = 500, steadystate_alt = upper);
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
@time beliefs, s = simulate_learning((options.burnin+1):options.N, s, beliefs, indices, options)

# Repeat learning with new shocks
noise_p = par.σ_p*randn((options.N - options.burnin))
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_p[(options.burnin+1):options.N] = simulate_ar(par.ρ_p, par.σ_p, options.N - options.burnin, noise_p)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs, s = simulate_learning((options.burnin+1):options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-5999):(500000-4999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.p[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))

"""
Plot phase diagram
"""
# Solution phase diagram
pyplot()
plot_ss(plot_points)
plot!(legend = :topright)
initial_ss = deepcopy(central)
starts = [(p=0.0,y=2.8,periods=16,arrows=[2]),
	(p=0.0,y=-3.8,periods=8,arrows=[2]),
	(p=0.0,y=-0.99,periods=100,arrows=[20]),
	(p=0.0,y=-0.9,periods=100,arrows=[45]),
	#(p=5.0,y=0.0,periods=100,arrows=[25]),
	#(p=-5.0,y=-0.0,periods=100,arrows=[25]),
	]
for start in starts
	initial_ss[:p] = start[:p]; initial_ss[:y] = start[:y];
	paths = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = start[:periods],
		magnitude = 1.0, persistence = par.ρ_y, show_plot = false)
	paths.y_lag = cat(start[:y],paths.y[1:(start.periods -1 )],dims=1)
	if start == starts[1]
		phase_arrow_plot(paths, [:y_lag,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 2:(start[:periods]), label = "Equilibrium paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths, [:y_lag,:p], arrow_points=start[:arrows], h_points = 2:(start[:periods]),
			v_points = 2:start[:periods], arrow_size = .5, final_arrow = true)
		end
end
plot!(size = (600,400))
savefig("figures/general/phase_compsource_sim.pdf")





















"""
Real source
"""

@everywhere par = (ϕ_pp = 0.95, ϕ_py = 0.5,
	ϕ_yy = 2.0, ϕ_yp = 0.25,
    α_2 = 0.0, α_3 = 0.1125,
	ρ_y = 0.5, σ_y = 0.5, ρ_p = 0.5, σ_p = 1.0);

a_11 = (1 - par.ϕ_py*par.ϕ_yp)/par.ϕ_pp
a_12 = - (par.ϕ_py*par.ϕ_yy)/par.ϕ_pp
a_21 = par.ϕ_yp
a_22 = par.ϕ_yy
A_matrix = [a_11 a_12;
	a_21 a_22]
eigen(A_matrix)

sstates = nlsolve(steady_states, [2.0, 2.0])
upper[:p] = sstates.zero[1];
upper[:y] = sstates.zero[2];
upper[:Ep_lead] = upper[:p];
# lower steady state
lower = Dict{Symbol,Float64}();
lower[:p] = sstates.zero[1];
lower[:y] = sstates.zero[2];
lower[:Ep_lead] = lower[:p];


"""
Perfect foresight phase
"""
pyplot()
plot_points = -4.0:0.01:4.0;
function plot_ss(plot_points)
	ss_plot = plot(xlabel = L"y_{t-1}", xlims = (minimum(plot_points),maximum(plot_points)),
    	ylabel = L"p_t", ylims = (-16.0,16.0),legend=:bottomright, yguidefontrotation=-90)
		plot!(plot_points,p_condition.(plot_points), label = "p condition", color = :black)
		display(plot!(y_condition.(5.0.*plot_points),5.0.*plot_points, label = "y condition", color = :green))
		# Plot perfect foresight paths
	return ss_plot
end

starts = [(p=0.05,y=-0.05,periods=20,arrows=[4]),
	(p=-0.05,y=0.05,periods=20,arrows=[4]),
	(p=0.5,y=-0.215128,periods=30,arrows=[12, 17]),
	(p=0.5,y=-0.215124,periods=30,arrows=[12, 20]),
	(p=-0.5,y=0.215128,periods=30,arrows=[12, 17]),
	(p=-0.5,y=0.215124,periods=30,arrows=[12, 20]),
	(p=0.5,y=-0.21512511575,periods=25,arrows=[12]),
	(p=-0.5,y=0.21512511575,periods=25,arrows=[12]),
	#(p=4.0,y=-1.5,periods=20,arrows=[4]),
	#(p=-4.0,y=1.5,periods=20,arrows=[4])
	]
plot_ss(plot_points)
plot!(legend = :topright)
initial_ss = deepcopy(central)
for start in starts
	initial_ss[:p] = start[:p]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods])

	if start == starts[1]
		phase_arrow_plot(paths1, [:y,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 1:(start[:periods]-1), label = "Perfect foresight paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths1, [:y,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 1:(start[:periods]-1), label = "", arrow_size = .5,
			final_arrow = true)
	end
end
plot!(size = (600,400))
savefig("figures/general/phase_realsource_perf.pdf")


"""
Re-train beliefs
"""
# Initialise beliefs
@everywhere beliefs = initialise_beliefs(options)
s = initialise_df(s, lower, gap = 500, steadystate_alt = upper)
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss)

# Simulate the learning for a set number of periods
noise_p = par.σ_p*randn(nrow(s))
s.ϵ_p = simulate_ar(par.ρ_p, par.σ_p, options.N, noise_p)
noise_y = par.σ_y*randn(nrow(s))
s.ϵ_y = simulate_ar(par.ρ_y, par.σ_y, options.N, noise_y)
plot(s.ϵ_y[1:1000])
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], lower, gap = 500, steadystate_alt = upper);
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
@time beliefs, s = simulate_learning((options.burnin+1):options.N, s, beliefs, indices, options)

# Repeat learning with new shocks
noise_p = par.σ_p*randn((options.N - options.burnin))
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_p[(options.burnin+1):options.N] = simulate_ar(par.ρ_p, par.σ_p, options.N - options.burnin, noise_p)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs, s = simulate_learning((options.burnin+1):options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-9999):(500000-8999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.p[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))

"""
Plot phase diagram
"""
# Solution phase diagram
pyplot()
plot_ss(plot_points)
plot!(legend = :topright)
initial_ss = deepcopy(central)
starts = [(p=0.0,y=-4.0,periods=20,arrows=[4]),
	(p=0.0,y=4.0,periods=20,arrows=[4]),
	(p=0.0,y=-0.65,periods=80,arrows=[2, 20]),
	(p=0.0,y=-0.79,periods=40,arrows=[2, 14]),
	#(p=5.0,y=0.0,periods=100,arrows=[25]),
	#(p=-5.0,y=-0.0,periods=100,arrows=[25]),
	]
for start in starts
	initial_ss[:p] = start[:p]; initial_ss[:y] = start[:y];
	paths = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = start[:periods],
		magnitude = 1.0, persistence = par.ρ_y, show_plot = false)
	paths.y_lag = cat(start[:y],paths.y[1:(start.periods -1 )],dims=1)
	if start == starts[1]
		phase_arrow_plot(paths, [:y_lag,:p], arrow_points=start[:arrows], h_points = 2:start[:periods],
			v_points = 2:(start[:periods]), label = "Equilibrium paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths, [:y_lag,:p], arrow_points=start[:arrows], h_points = 2:(start[:periods]),
			v_points = 2:start[:periods], arrow_size = .5, final_arrow = true)
		end
end
plot!(size = (600,400))
savefig("figures/general/phase_realsource_sim.pdf")







































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
savefig("figures/general/irf1_general_sim.pdf")
plot(paths2.π, label =L"\pi_t", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths2.y, label =L"y_t");plot!(paths2.ϵ_y, label =L"\epsilon_{y,t}")
plot!(size=(300,200))
savefig("figures/general/irf2_general_sim.pdf")
plot(paths3.π, label =L"\pi_t", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths3.y, label =L"y_t");plot!(paths3.ϵ_y, label =L"\epsilon_{y,t}")
plot!(size=(300,200))
savefig("figures/general/irf3_general_sim.pdf")
plot(paths4.π, label =L"\pi_t", xlabel = "Periods", legend = false,ylims = (-3.0,3.0))
plot!(paths4.y, label =L"y_t");plot!(paths4.ϵ_y, label =L"\epsilon_{y,t}")
plot!(size=(300,200))
savefig("figures/general/irf4_general_sim.pdf")

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
savefig("figures/general/irfs_general_sim.pdf")



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
savefig("figures/general/heatmap_general_sim.pdf")


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
savefig("figures/general/heatmap_errors_general_sim.pdf")



Rsq_π = 1 - var(heatmap_df.π-heatmap_df.Eπ)/var(heatmap_df.π)
Rsq_y = 1 - var(heatmap_df.y-heatmap_df.Ey)/var(heatmap_df.y)

















































"""
End of script
"""
