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

@everywhere options = EcoNNetOptions(infoset = [:π_lag, :y_lag, :ϵ_π, :ϵ_y],
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

@everywhere par = (β = 0.95, κ = 0.05, η = 0.95, σ = 0.25,
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
Calculate steady states
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

#for tt in 1:50
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=-1.8,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=1.8,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-2.2,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=2.2,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-2.5,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=2.5,y=3.5,periods=100,arrows=[10,50,98])
	]
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
savefig("figures/pw_linear/perf_phase_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")


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
options.burnin = 100000;
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], lower, gap = 500, steadystate_alt = upper);
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
options.burnin_use_net = false;
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
@time beliefs,s = simulate_learning(options.burnin:options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-5999):(500000-4999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.π[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))
savefig("figures/pw_linear/sim_series_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")


export_df = s[options.N-99999:options.N,:]
rename!(export_df, replace.(names(export_df), "π" => "pi"))
rename!(export_df, replace.(names(export_df), "ϵ" => "epsilon"))
export_df.r = Taylor_condition.(export_df.pi)
export_df = export_df[:,[:epsilon_pi, :epsilon_y, :pi, :y, :r]]

CSV.write("estimation/pwlin/pwlin_sim_pistar"*rep_pnt(par.π_star)*".csv",
	export_df)

"""
Plot phase diagram
"""

pyplot()
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=3.0,y=3.0,periods=100,arrows=[10]),
	(π=-3.,y=-3.,periods=100,arrows=[10]),
	(π=0.1,y=0.1,periods=100,arrows=[22, 65]),
	(π=-0.1,y=-0.1,periods=100,arrows=[22, 65]),
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
savefig("figures/pw_linear/phase_nnet_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")








"""
Smaller π_star (0.5)
"""

@everywhere par = (β = 0.95, κ = 0.05, η = 0.95, σ = 0.25,
	ϕ_π = 0.5, π_star = 0.5, α = 0.75,
	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);

"""
Perfect foresight phase
"""
plot_points = -4.0:0.01:4.0;
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=-1.5,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=1.5,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-2.0,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=2.0,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-1.85,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=1.85,y=3.5,periods=100,arrows=[10,50,98])
	]
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
savefig("figures/pw_linear/perf_phase_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")


"""
Re-train beliefs
"""
@everywhere beliefs = initialise_beliefs(options)

# Simulate the learning for a set number of periods
noise_π = par.σ_π*randn((options.N - options.burnin))
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_π[(options.burnin+1):options.N] = simulate_ar(par.ρ_π, par.σ_π, options.N - options.burnin, noise_π)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs,s = simulate_learning(options.burnin:options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-6999):(500000-4999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.π[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))
savefig("figures/pw_linear/sim_series_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")



export_df = s[options.N-99999:options.N,:]
rename!(export_df, replace.(names(export_df), "π" => "pi"))
rename!(export_df, replace.(names(export_df), "ϵ" => "epsilon"))
export_df.r = Taylor_condition.(export_df.pi)
export_df = export_df[:,[:epsilon_pi, :epsilon_y, :pi, :y, :r]]

CSV.write("estimation/pwlin/pwlin_sim_pistar"*rep_pnt(par.π_star)*".csv",
	export_df)


"""
Plot phase diagram
"""

pyplot()
plot_ss(plot_points)
initial_ss = deepcopy(central)
periods = 200
starts = [(π=3.0,y=3.0,periods=periods,arrows=[10,periods-100]),
	(π=-3.,y=-3.,periods=periods,arrows=[10,periods-100]),
	(π=0.1,y=0.1,periods=periods,arrows=[102,periods-2]),
	(π=-0.1,y=-0.1,periods=periods,arrows=[102,periods-2]),
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
savefig("figures/pw_linear/phase_nnet_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")





















"""
Smallest π_star (0.25)
"""

@everywhere par = (β = 0.95, κ = 0.05, η = 0.95, σ = 0.25,
	ϕ_π = 0.5, π_star = 0.25, α = 0.75,
	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);

"""
Perfect foresight phase
"""
plot_points = -4.0:0.01:4.0;
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=-1.6,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=1.6,y=3.5,periods=100,arrows=[10,50,98]),
	(π=-1.62,y=-3.5,periods=200,arrows=[10,50,198]),
	(π=1.62,y=3.5,periods=200,arrows=[10,50,198]),
	(π=-1.7,y=-3.5,periods=100,arrows=[10,50,98]),
	(π=1.7,y=3.5,periods=100,arrows=[10,50,98])
	]
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
savefig("figures/pw_linear/perf_phase_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")


"""
Re-train beliefs
"""
@everywhere beliefs = initialise_beliefs(options)

# Simulate the learning for a set number of periods
noise_π = par.σ_π*randn((options.N - options.burnin))
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_π[(options.burnin+1):options.N] = simulate_ar(par.ρ_π, par.σ_π, options.N - options.burnin, noise_π)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs,s = simulate_learning(options.burnin:options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-6999):(500000-4999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.π[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))
savefig("figures/pw_linear/sim_series_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")

"""
Plot phase diagram
"""

pyplot()
plot_ss(plot_points)
initial_ss = deepcopy(central)
periods = 5000
starts = [(π=3.0,y=3.0,periods=periods,arrows=[10,periods-4900]),
	(π=-3.,y=-3.,periods=periods,arrows=[10,periods-4900]),
	(π=0.2,y=0.2,periods=periods,arrows=[102,periods-2500]),
	(π=-0.2,y=-0.2,periods=periods,arrows=[102,periods-4600]),
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
savefig("figures/pw_linear/phase_nnet_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")

















"""
Larger π_star 2.0
"""

@everywhere par = (β = 0.95, κ = 0.05, η = 0.95, σ = 0.25,
	ϕ_π = 0.5, π_star = 2.0, α = 0.75,
	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.5, σ_π = 0.2);

"""
Perfect foresight phase
"""
plot_points = -6.0:0.01:6.0;
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=-4.25,y=-5.5,periods=100,arrows=[10,50,98]),
	(π=4.25,y=5.5,periods=100,arrows=[10,50,98]),
	(π=-2.6,y=-5.5,periods=200,arrows=[10,50,198]),
	(π=2.6,y=5.5,periods=200,arrows=[10,50,198]),
	(π=-1.7,y=-5.5,periods=100,arrows=[10,50,98]),
	(π=1.7,y=5.5,periods=100,arrows=[10,50,98])
	]
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
savefig("figures/pw_linear/perf_phase_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")


"""
Re-train beliefs
"""
@everywhere beliefs = initialise_beliefs(options)

# Simulate the learning for a set number of periods
noise_π = par.σ_π*randn((options.N - options.burnin))
noise_y = par.σ_y*randn((options.N - options.burnin))
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_π[(options.burnin+1):options.N] = simulate_ar(par.ρ_π, par.σ_π, options.N - options.burnin, noise_π)
s.ϵ_y[(options.burnin+1):options.N] = simulate_ar(par.ρ_y, par.σ_y, options.N - options.burnin, noise_y)
@time beliefs,s = simulate_learning(options.burnin:options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (500000-6999):(500000-4999)
plot(layout=(2,1),legend = false,  link = :x)
plot!(s.π[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))
savefig("figures/pw_linear/sim_series_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")

"""
Plot phase diagram
"""

pyplot()
plot_ss(plot_points)
initial_ss = deepcopy(central)
periods = 5000
starts = [(π=5.0,y=5.0,periods=periods,arrows=[10,periods-4900]),
	(π=-5.,y=-5.,periods=periods,arrows=[10,35,periods-4900]),
	(π=0.2,y=0.2,periods=periods,arrows=[102,periods-2500]),
	(π=-0.2,y=-0.2,periods=periods,arrows=[102,periods-4600])
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
savefig("figures/pw_linear/phase_nnet_pistar"*rep_pnt(par.π_star)*
	"_alpha"*rep_pnt(par.α)*".pdf")
