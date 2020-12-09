cd("/Users/julianashwin/Documents/GitHub/EcoNNet.jl")

"""
Add extra cores if you want to parallelise later
"""

using LaTeXStrings, TableView
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

save_figs = false
create_gifs = false

@everywhere options = EcoNNetOptions(infoset = [:π_lag, :y_lag, :ϵ_π, :ϵ_y],
    expectations = [:Eπ,:Ey,:Eπ_lead],
    endogenous = [:π, :y],
    exogenous = [:ϵ_π,:ϵ_y],
    states = [:y_lag, :ϵ_π, :ϵ_y],
	auxiliary = [:r],
    N = 50000, num_nodes = 8, activation = relu, window = 400)

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
central[:π] = sstates.zero[1]; central[:Eπ] = central[:π];
central[:y] = sstates.zero[2]; central[:Ey] = central[:y];
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
	ss_plot = plot(xlabel = L"y_{t-1}", xlims = (-4.0,4.0),
    	ylabel = L"\pi_t", ylims = (-4.0, 4.0),legend=:bottomright, yguidefontrotation=-90)
		plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
		display(plot!(IS_condition.(plot_points),plot_points, label = "IS Curve", color = :green))
		# Plot perfect foresight paths
	return ss_plot
end

#for tt in 1:50
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=0.75,y=-1.5,periods=100,arrows=[10,60]),
	(π=-0.75,y=1.5,periods=100,arrows=[10,60]),
	(π=0.0,y=-1.5,periods=62,arrows=[10,60]),
	(π=0.0,y=1.5,periods=62,arrows=[10,60]),
	(π=-1.024,y=-3.0,periods=100,arrows=[10,55]),
	(π=-2.350785,y=-3.0,periods=70,arrows=[25]),
	(π=1.024,y=3.0,periods=100,arrows=[10,55]),
	(π=2.350785,y=3.0,periods=70,arrows=[25])]
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
if save_figs
	plot!(size = (600,400))
	savefig("figures/phase_illus_perf.pdf")
end

if create_gifs
starts = [(π=0.75,y=-1.5,periods=100,arrows=[10,60]),
	(π=-0.75,y=1.5,periods=100,arrows=[10,60]),
	(π=0.0,y=-1.5,periods=100,arrows=[10,60]),
	(π=0.0,y=1.5,periods=100,arrows=[10,60]),
	(π=-1.024,y=-3.0,periods=80,arrows=[10,55]),
	(π=-2.350785,y=-3.0,periods=60,arrows=[25]),
	(π=1.024,y=3.0,periods=80,arrows=[10,55]),
	(π=2.350785,y=3.0,periods=60,arrows=[25])]

periods = 100

anim = @animate for tt in 2:periods
	plot_ss(plot_points)
	plot!(size = (600,400))
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
gif(anim, "figures/phase_illus_perf.gif", fps = 15)
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
s = initialise_df(s, ss);
s = initialise_df(s, central);
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss,cutoff= false)


"""
Run learning simulation
"""
noise_π = par.σ_π*randn(nrow(s))
s.ϵ_π = simulate_ar(par.ρ_π, par.σ_π, options.N, noise_π)
noise_y = par.σ_y*randn(nrow(s))
s.ϵ_y = simulate_ar(par.ρ_y, par.σ_y, options.N, noise_y)
plot(s.ϵ_y[1:200])
options.burnin = 1000;
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
options.burnin_use_net = false;
options.learning_gap = 1;
options.plotting_gap = 1;
options.plot_vars = [:π, :y, :Eπ, :Ey]

# Simulate the learning for a set number of periods
gr() # Set GR backend for plots as it's the fastest
@time beliefs,s = simulate_learning(options.burnin:1900, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
plot_range = (1900-999):(1900)
plot(layout=(2,1),legend = :topright,  link = :x)
plot!(s.π[plot_range], subplot = 1, label = L"\pi_t", ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.Eπ[plot_range], subplot = 1, label = L"E \pi_t", linestyle = :dash, yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, label = L"y_t", ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(s.Ey[plot_range], subplot = 2, label = L"Ey_t", linestyle = :dash, yguidefontrotation=-90)
plot!(size = (600,300))
savefig("figures/illus_realtime_series.pdf")


"""
Plot gif of real time learning
"""

start_point = 901
periods = 900
anim = @animate for tt in 1:periods
	plot_range = start_point:(start_point + 100+  tt)
	plot(layout=(2,1),legend = :topright,  link = :x, xlims = (1, 1000))
	plot!(s.π[plot_range], subplot = 1, label = L"\pi_t", ylabel = L"\pi_t",
		ylims = (-3.3,5), yguidefontrotation=-90)
	plot!(s.Eπ[plot_range], subplot = 1, label = L"E \pi_t", linestyle = :dash, yguidefontrotation=-90)
	plot!(s.y[plot_range], subplot = 2, label = L"y_t", ylabel = L"y_t",
		ylims = (-4, 3), yguidefontrotation=-90, xlabel = "Periods")
	plot!(s.Ey[plot_range], subplot = 2, label = L"Ey_t", linestyle = :dash, yguidefontrotation=-90)
	display(plot!(size = (600,300)))

end
gif(anim, "figures/illus_realtime_series.gif", fps = 100)



"""
End of script
"""
