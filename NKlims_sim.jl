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

@everywhere options = EcoNNetOptions(infoset = [:π_lag, :y_lag, :r_lag, :ϵ_π, :ϵ_y],
    expectations = [:Eπ,:Ey,:Eπ_lead, :Ey_lead],
    endogenous = [:π, :y, :r],
    exogenous = [:ϵ_π,:ϵ_y],
    states = [:π_lag, :y_lag, :r_lag, :ϵ_π, :ϵ_y],
    N = 500000, num_nodes = 32, activation = relu, window = 40000)

@everywhere beliefs = initialise_beliefs(options)

"""
Define the parameters as a Named Tuple.
"""

@everywhere par = (β = 0.95, κ = 0.1,
	η = 0.5, σ = 0.25, ϕ_π = 1.5, ϕ_y = 0.01,
	ρ_π = 0.5, σ_π = 0.001, ρ_y = 0.5, σ_y = 0.001,
	R_lim = -0.02, π_lim = -0.03, y_lim = -0.05);

"""
State the equilibrium conditions of the model as a function which returns
    a vector of zeros
"""

@everywhere function equilibrium_conditions_fast(F::Array{Float64,1},
    x::Array{Float64,1},states_input::Array{Float64,1},predictions_input::Array{Float64,1})
    # Manually unpack the states
    π_lag::Float64 = states_input[1]
    y_lag::Float64 = states_input[2]
	r_lag::Float64 = states_input[3]
    ϵ_π::Float64 = states_input[4]
	ϵ_y::Float64 = states_input[5]
    # and the predictions
    #Ep::Float64 = predictions_input[1]
    #Ey::Float64 = predictions_input[2]
    Eπ_lead::Float64 = predictions_input[3]
	Ey_lead::Float64 = predictions_input[4]
    # and the endogenous variables
    π::Float64 = x[1]
    y::Float64 = x[2]
	r::Float64 = x[3]

    # p[t] = ϕ_pp*pe[t+1] + ϕ_py*y[t] + α_2*y[t]^2 - α_3*y[t]^3 + ϵ_p[t]
    F[1] =  max(par.β*Eπ_lead + par.κ*y, par.π_lim)  + ϵ_π - π ::Float64
    # y[t] = η*y[t-1] + ϕ_yp*p[t] + ϵ_y[t]
	F[3] = r - max(par.ϕ_π*π + par.ϕ_y*y, par.R_lim)
    F[2] = max((1-par.η)*Ey_lead + par.η*y_lag - par.σ*(r - Eπ_lead), par.y_lim)+ ϵ_y - y ::Float64

    return F
end


"""
Define equations of motion under perfect foresight as a useful comparison
	Input - first lagged state variables, then endogenous variables
"""
function perfect_foresight(inputs)
    # Manually unpack the states (endogenous variables in same order as in options)
	# Lagged state variables first, then current endogenous
    π_lag = max(inputs[1], par.π_lim)
    y_lag = max(inputs[2], par.y_lim)
	r_lag = max(inputs[3], par.R_lim)
	π_t = max(inputs[4], par.π_lim)
	y_t = max(inputs[5], par.y_lim)
	r_t = max(inputs[6], par.R_lim)



    # π[t+1] = 1/β*(π[t] - κ*y[t])
    π_new = max((1/par.β)*(π_t - par.κ*y_t),par.π_lim)
    # y[t+1] = η*y[t] - σ*(ϕ_π*π[t+1] + α π[t+1]^3 - π[t+2])
	r_t = max(par.ϕ_π*π_t + par.ϕ_y*y_t, par.R_lim)
	y_new = (1/(1-par.η))*(y_t - par.η*y_lag + par.σ*(
		r_t - π_new))
	y_new = max(y_new,par.y_lim)
	r_new = max(par.ϕ_π*π_new + par.ϕ_y*y_new, par.R_lim)
	# Impose upper and lower bounds to allow plotting
	π_new = min(π_new,max(π_new,-1e6),1e6)
	y_new = min(y_new,max(y_new,-1e6),1e6)
	r_new = min(r_t,max(r_t,-1e6),1e6)

    outputs = [π_new, y_new, r_new]

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
	r::Float64 = x[3]
    F[1]::Float64 = max(par.β*π + par.κ*y,par.π_lim) - π
    F[2]::Float64 = max((1-par.η)*y + par.η*y - par.σ*(
		max(par.ϕ_π*π + par.ϕ_y*y, par.R_lim) - π),par.y_lim) - y
	F[3]::Float64 = max(par.ϕ_π*π + par.ϕ_y*y, par.R_lim) - r
    return F
end

# Exogenous processes
ss = Dict{Symbol,Float64}()
ss[:ϵ_π] = 0.0
ss[:ϵ_y] = 0.0
# upper steady state
sstates = nlsolve(steady_states, [2.0, 2.0, 2.0])
upper = Dict{Symbol,Float64}();
upper[:π] = sstates.zero[1];
upper[:y] = sstates.zero[2];
upper[:r] = sstates.zero[3];
upper[:Eπ_lead] = upper[:π];
# central steady state
sstates = nlsolve(steady_states, [-0.02, -0.02, -0.02])
central = Dict{Symbol,Float64}();
central[:π] = sstates.zero[1];
central[:y] = sstates.zero[2];
central[:r] = sstates.zero[3];
central[:Eπ_lead] = central[:π];
# lower steady state
lower = Dict{Symbol,Float64}();
sstates = nlsolve(steady_states, [-2.0, -2.0, -2.0])
lower[:π] = sstates.zero[1];
lower[:y] = sstates.zero[2];
lower[:r] = sstates.zero[3];
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
pyplot()
plot_points = -1:0.001:1
function plot_ss(plot_points)
	ss_plot = plot(xlabel = L"y_t", xlims = (-0.07,0.07),
	    ylabel = L"\pi_t", ylims = (-0.07, 0.07),legend=:topright)
	plot!(plot_points, ZLB_bind_condition.(plot_points),
		fillrange=[minimum(plot_points)*ones(len(plot_points))], fillalpha = 0.5,
		 color = :paleturquoise1, label = "ELB")
	plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
	plot!(par.y_lim:0.001:EE_kink,EE_normal_condition.(par.y_lim:0.001:EE_kink),
	 	label = "Euler Equation", color = :green)
	plot!(par.y_lim:0.001:EE_kink,EE_ZLB_condition.(par.y_lim:0.001:EE_kink),
	 	label = "", color = :green)
	plot!([par.y_lim,par.y_lim],[par.R_lim,-0.07],
	 	label = "", color = :green)
	return ss_plot
end
## Plot perfect foresight paths
pyplot()
phase_plot = plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=0.003,y=0.003, periods=5,arrows=[4]),
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
	print(start)
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = pf_path(initial_ss, periods = start[:periods], lags = 2)
	if start == starts[1]
		phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
			v_points = 1:start[:periods], label = "Perfect foresight", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
			v_points = 1:start[:periods], arrow_size = .5,
			final_arrow = true)
	end
end
display(phase_plot)
plot!(size = (600,400))

savefig("figures/NK_model/phase_NKlim_perf.pdf")

for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = 100,
		magnitude = start[:ϵ_y], persistence = par.ρ_y, show_plot = false)


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
noise_y = par.σ_y*randn(options.N) + alternating_shocks(options.N, gap=200, mag = 0.02)
s.ϵ_y = simulate_ar(par.ρ_y, 0.004, options.N, noise_y)
s.ϵ_π = simulate_ar(par.ρ_π, 0.004, options.N, Normal())
plot(s.ϵ_y[1:2000])
plot(s.ϵ_π[1:2000])
options.burnin = 50000;
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], lower, gap = 500, steadystate_alt = upper);
s[1:options.burnin,:] = initialise_df(s[1:options.burnin,:], ss);
options.burnin_use_net = false;
options.learning_gap = 50000;
options.plotting_gap = 50000;
options.window = 49999;
options.plot_vars = [:π, :y, :r, :Eπ, :Ey]

# Simulate the learning for a set number of periods
gr() # Set GR backend for plots as it's the fastest
s[1:options.burnin,:]= s[(options.N - options.burnin + 1):options.N,:]
@time beliefs,s = simulate_learning(options.burnin:options.N, s, beliefs, indices, options)

# Plot simulated time series
gr()
plot_range = (90000):(91000)
plot(layout=(2,1), xlabel = "Periods",legend = false)
plot!(s.π[plot_range], subplot = 1, ylabel = "Inflation")
plot!(s.y[plot_range], subplot = 2, ylabel = "Output")
plot!(size = (1000,400), left_margin = 7mm)
savefig("figures/NK_sim_series.png")



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
plot_points = -1:0.001:1
function plot_ss(plot_points)
	ss_plot = plot(xlabel = L"y_t", xlims = (-0.07,0.07),
	    ylabel = L"\pi_t", ylims = (-0.07, 0.07),legend=:topright)
	plot!(plot_points, ZLB_bind_condition.(plot_points),
		fillrange=[minimum(plot_points)*ones(len(plot_points))], fillalpha = 0.5,
		 color = :paleturquoise1, label = "ELB")
	plot!(plot_points,NKPC_condition.(plot_points), label = "Phillips Curve", color = :black)
	plot!(par.y_lim:0.001:EE_kink,EE_normal_condition.(par.y_lim:0.001:EE_kink),
	 	label = "Euler Equation", color = :green)
	plot!(par.y_lim:0.001:EE_kink,EE_ZLB_condition.(par.y_lim:0.001:EE_kink),
	 	label = "", color = :green)
	plot!([par.y_lim,par.y_lim],[par.R_lim,-0.07],
	 	label = "", color = :green)
	return ss_plot
end

## Plot equilibrium phase diagram
phase_plot = plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=-0.04,y=-0.06,ϵ_y=-0.0,periods=2,arrows=[1]),
	(π=-0.04,y=-0.027,ϵ_y=0.0,periods=15,arrows=[7]),
	(π=-0.04,y=-0.029,ϵ_y=0.0,periods=24,arrows=[10]),
	(π=0.0,y=-0.032,ϵ_y=0.0,periods=40,arrows=[]),
	(π=0.0,y=-0.034,ϵ_y=0.0,periods=60,arrows=[6]),
	(π=-0.0175,y=0.035,ϵ_y=0.0,periods=12,arrows=[]),
	(π=0.03,y=0.035,ϵ_y=0.0,periods=12,arrows=[3]),
	]

for start in starts
	initial_ss[:π] = start[:π]; initial_ss[:y] = start[:y];
	paths1 = irf(:ϵ_y, initial_ss, beliefs, shock_period = 2, periods = 100,
		magnitude = start[:ϵ_y], persistence = par.ρ_y, show_plot = false)
	if start == starts[1]
		phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
			v_points = 1:start[:periods], label = "Equilibrium paths", arrow_size = .5,
			final_arrow = true)
	else
		phase_arrow_plot(paths1, [:y,:π], arrow_points=start[:arrows], h_points = 1:start[:periods],
			v_points = 1:start[:periods], arrow_size = .5,
			final_arrow = true)
	end

end
display(phase_plot)
plot!(size = (600,400))
savefig("figures/NK_model/NK_phase_sim.pdf")











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
Compare PLM and ALM
"""
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
