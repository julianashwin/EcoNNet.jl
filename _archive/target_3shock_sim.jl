cd("/Users/julianashwin/Documents/GitHub/EcoNNet.jl")

"""
Add extra cores if you want to parallelise later
"""

using LaTeXStrings, TableView, CSV, Setfield
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

@everywhere options = EcoNNetOptions(
	infoset = [:ϵ_π],
    expectations = [:Eπ,:Ey,:Eπ_lead],
    endogenous = [:π, :y],
    exogenous = [:ϵ_π,:ϵ_y],
    states = [:ϵ_π, :ϵ_y],
	auxiliary = [:r],
    N = 50000, num_nodes = 2, activation = relu, window = 500,
	burnin = 5000, learning_gap = 1000, plotting_gap = 100)

@everywhere beliefs = initialise_beliefs(options)
params(beliefs)

"""
Define the parameters as a Named Tuple.
"""

@everywhere par = (β = 1, κ = 1, η = 0, σ = 1,
	ϕ_π = 0, π_star = 2.0, α = 1.5, ϵ = 5.0,
	ρ_y = 0.5, σ_y = 0.2, ρ_π = 0.0, σ_π = 0.0);

"""
State the equilibrium conditions of the model as a function which returns
    a vector of zeros
"""

@everywhere function equilibrium_conditions_fast(F::Array{Float64,1},
    x::Array{Float64,1},states_input::Array{Float64,1},predictions_input::Array{Float64,1})
    # Manually unpack the states
    #p_lag::Float64 = states_input[1]
    y_lag::Float64 = 0.#states_input[1]
    ϵ_π::Float64 = states_input[1]
	ϵ_y::Float64 = states_input[2]
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
How are the equilibrium conditions actually solved?
"""
ys = -6:0.1:6
πs = -6:0.1:6

# IS condition (including the Taylor Rule)
function IS_curve(π_, Eπ, ϵ)
	if (abs(π_) <= par.π_star)
		y = Eπ - 0
	elseif π_ > par.π_star
		y = Eπ - par.α * (π_ - par.π_star)
	elseif π_ < par.π_star
		y = Eπ - par.α * (π_ + par.π_star)
	end
	return y
end
# NKPC condition (including the Taylor Rule)
function NKPC_curve(y, Eπ, ϵ)
	π_ = Eπ + y + ϵ
	return π_
end


pyplot()
plot(xlabel = L"y_{t}", ylabel = L"\pi_t", legend=:topleft, yguidefontrotation=-90,
	xlims = (-6,6), ylims = (-6,6))
# Plot IS condition
plot!(IS_curve.(πs, [0.], [0.,]), πs, label = "IS", color = :green, linestyle = :dot)
# Plot NKPC condition
plot!(ys,NKPC_curve.(ys, [0.], [5.]), label = "NKPC", color = :blue,
	annotations = ([2.2], [5.], L"E\pi_{t+1} = 0, \epsilon_t = 5", :blue))
plot!(ys,NKPC_curve.(ys, [0.], [0]), label = "", color = :blue,
	annotations = ([2.2], [0], L"E\pi_{t+1} = 0, \epsilon_t = 0", :blue))
plot!(ys,NKPC_curve.(ys, [0.], [-5.0]), label = "", color = :blue,
	annotations = ([2.4], [-5.], L"E\pi_{t+1} = 0, \epsilon_t = -5", :blue))
plot!(size = (600,400))
savefig("figures/analytics/eqlm_conditions_example.pdf")



function single_condition(π_, Eπ, ϵ)
	if (abs(π_) <= par.π_star)
		y = Eπ - 0
	elseif π_ > par.π_star
		y = Eπ - par.α * (π_ - par.π_star)
	elseif π_ < par.π_star
		y = Eπ - par.α * (π_ + par.π_star)
	end
	π_ = Eπ + y + ϵ
	return π_
end

pyplot()
plot(xlabel = L"π_{t}", ylabel = L"\pi_t", legend=:topleft, yguidefontrotation=-90,
	xlims = (-6,6), ylims = (-6,6))
# Plot IS condition
plot!(πs, πs, label = "45 degree", color = :green, linestyle = :dot)
# Plot NKPC condition
plot!(ys,single_condition.(πs, [0.], [5.0]), label = "Eqlm condition", color = :blue,
	annotations = ([0.], [5.5], L"E\pi_{t+1} = 0, \epsilon_t = 5", :blue))
plot!(ys,single_condition.(ys, [0.], [0.]), label = "", color = :blue,
	annotations = ([0.], [0.5], L"E\pi_{t+1} = 0, \epsilon_t = 0", :blue))
plot!(ys,single_condition.(ys, [0.], [-5]), label = "", color = :blue,
	annotations = ([0.], [-4.5], L"E\pi_{t+1} = 0, \epsilon_t = -5", :blue))
plot!(size = (600,400))
savefig("figures/analytics/single_eqlm_condition_example.pdf")



"
Characterise the REE
"
P = [0.5 0.25 0.25; 0.25 0.5 0.25; 0.25 0.25 0.5]
P2 = [0.6 0.3 0.1; 0.2 0.6 0.2; 0.1 0.3 0.6]
function L_condition_U(L, M, P)
	if (abs(L) <= par.π_star)
		U = (-1/(2*P[1,3]))*(2*(P[1,1]* L + P[1,2]*M ) - L - par.ϵ)
	elseif L > par.π_star
		U = (-1/(2*P[1,3]))*(2*(P[1,1]* L + P[1,2]*M ) - par.α*(L - par.π_star) - L - par.ϵ)
	elseif L < par.π_star
		U = (-1/(2*P[1,3]))*(2*(P[1,1]* L + P[1,2]*M ) - par.α*(L + par.π_star) - L - par.ϵ)
	end
	return U
end
function U_condition_L(M, U, P)
	if (abs(U) <= par.π_star)
		L = (-1/(2*P[3,1]))*(2*(P[3,3]* U + P[3,2]*M ) - U + par.ϵ)
	elseif U > par.π_star
		L = (-1/(2*P[3,1]))*(2*(P[3,3]* U + P[3,2]*M ) - par.α*(U - par.π_star) - U + par.ϵ)
	elseif U < par.π_star
		L = (-1/(2*P[3,1]))*(2*(P[3,3]* U + P[3,2]*M ) - par.α*(U + par.π_star) - U + par.ϵ)
	end
	return L
end
function U_condition_M(L, U, P)
	if (abs(U) <= par.π_star)
		M = (-1/(2*P[3,2]))*(2*(P[3,3]* U + P[3,1]*L ) - U + par.ϵ)
	elseif U > par.π_star
		M = (-1/(2*P[3,2]))*(2*(P[3,3]* U + P[3,1]*L ) - par.α*(U - par.π_star) - U + par.ϵ)
	elseif U < par.π_star
		M = (-1/(2*P[3,2]))*(2*(P[3,3]* U + P[3,1]*L ) - par.α*(U + par.π_star) - U + par.ϵ)
	end
	return M
end
function M_condition_L(M, U, P)
	if (abs(M) <= par.π_star)
		L = (-1/(2*P[2,1]))*(2*(P[2,2]*M +  P[2,3]*U) - M)
	elseif M > par.π_star
		L = (-1/(2*P[2,1]))*(2*(P[2,2]*M +  P[2,3]*U) - par.α*(M - par.π_star) - M)
	elseif M < par.π_star
		L = (-1/(2*P[2,1]))*(2*(P[2,2]*M +  P[2,3]*U) - par.α*(M + par.π_star) - M)
	end
	return L
end


# If we impose that M = 0, we can graphically show the solutions
Us = -30:1:30
Ls = -30:1:30
Ms = -30:1:30
## With baseline P martrix
par = @set par.ϵ = 10.0
pyplot()
plot(xlabel = L"L = \pi_t | ϵ_t = -\epsilon", ylabel = L"U = \pi_t | ϵ_t = \epsilon",
	legend=:topleft ,xlims = (-30,30), ylims = (-30,30))
plot!(zeros(2), [-30, 30], color = :black, linestyle = :dot, label = "")
plot!([-30, 30], zeros(2), color = :black, linestyle = :dot, label = "")
# Plot L condition
plot!(Ls, L_condition_U.(Ls, [-20.], [P]), label = "M = -20", color = :blue)
plot!(Ls, L_condition_U.(Ls, [0.], [P]), label = "M = 0", color = :purple)
plot!(Ls, L_condition_U.(Ls, [20.], [P]), label = "M = +20", color = :red)
# Plot U condition
plot!(U_condition_L.([-20.], Us, [P]), Us, label = "", color = :blue)
plot!(U_condition_L.([0.], Us, [P]), Us, label = "", color = :purple)
plot!(U_condition_L.([20.], Us, [P]), Us, label = "", color = :red)
# Highlight the equilibria
scatter!([-18.], [-8.], color = :blue, label = "", markersize = 5)
scatter!([-6.5], [6.5], color = :purple, label = "", markersize = 5)
scatter!([8.], [18.], color = :red, label = "", markersize = 5)
plot!(size = (600,400))
savefig("figures/analytics/UL_givenM.pdf")



## With alternative P matrix
pyplot()
plot(xlabel = L"L = \pi_t | ϵ_t = -\epsilon", ylabel = L"U = \pi_t | ϵ_t = \epsilon",
	legend=:topleft ,xlims = (-30,30), ylims = (-30,30))
plot!(zeros(2), [-30, 30], color = :black, linestyle = :dot, label = "")
plot!([-30, 30], zeros(2), color = :black, linestyle = :dot, label = "")
# Plot L condition
plot!(Ls, L_condition_U.(Ls, [-20.], [P2]), label = "M = -20", color = :blue)
plot!(Ls, L_condition_U.(Ls, [0.], [P2]), label = "M = 0", color = :purple)
plot!(Ls, L_condition_U.(Ls, [20.], [P2]), label = "M = +20", color = :red)
# Plot U condition
plot!(U_condition_L.([-20.], Us, [P2]), Us, label = "", color = :blue)
plot!(U_condition_L.([0.], Us, [P2]), Us, label = "", color = :purple)
plot!(U_condition_L.([20.], Us, [P2]), Us, label = "", color = :red)



## Solve
@everywhere function analytical_REE(F::Array{Float64,1}, x::Array{Float64,1})

	P = [0.5 0.25 0.25; 0.25 0.5 0.25; 0.25 0.25 0.5]

	L::Float64 = x[1]
    M::Float64 = x[2]
	U::Float64 = x[3]
    F[1] = L - M_condition_L(M, U, P) ::Float64
	F[2] = M - U_condition_M(L, U, P) ::Float64
	F[3] = U - L_condition_U(L, M, P) ::Float64
    return F
end

nlsolve(analytical_REE, [0.0, 0.0, 0.0])


## Three dimensions
soln_vars = [:soln_no, :π_star, :ϵ, :ϵ_π, :y, :π, :r, :Eπ_lead]
REE_solns = DataFrame(zeros(3, length(soln_vars)), soln_vars);
# Solve for L


P = [0.6 0.3 0.1; 0.2 0.6 0.2; 0.1 0.3 0.6]
P = [0.5 0.25 0.25; 0.25 0.5 0.25; 0.25 0.25 0.5]
par = @set par.ϵ = 10.0 

Us = -10:0.1:10
Ls = -10:0.1:10
Ms = -10:0.1:10
state_labels = [(collect(x)) for x in (vec(collect(Iterators.product(Ls, Ms, Us))))]
state_grid = permutedims(reshape(hcat(state_labels...), (length(state_labels[1]), length(state_labels))))

comb_conditions = DataFrame(state_grid, [:L, :M, :U])
comb_conditions.L_soln = M_condition_L.(comb_conditions.M,comb_conditions.U,[P])
comb_conditions.M_soln = U_condition_M.(comb_conditions.L,comb_conditions.U,[P])
comb_conditions.U_soln = L_condition_U.(comb_conditions.L,comb_conditions.M,[P])

solns = comb_conditions[(abs.(comb_conditions.L .- comb_conditions.L_soln) .< 0.05) .*
						(abs.(comb_conditions.M .- comb_conditions.M_soln) .< 0.05) .*
						(abs.(comb_conditions.U .- comb_conditions.U_soln) .< 0.05) ,:]

						
plot(xlabel = L"\epsilon_t", ylabel = L"\pi_t")
plot!([-par.ϵ, 0, par.ϵ], collect(solns[1,[:L, :M, :U]]), color = 1, label = "")
scatter!([-par.ϵ, 0, par.ϵ], collect(solns[1,[:L, :M, :U]]), color = 1, label = "")
plot!([-par.ϵ, 0, par.ϵ], collect(solns[2,[:L, :M, :U]]), color = 2, label = "")
scatter!([-par.ϵ, 0, par.ϵ], collect(solns[2,[:L, :M, :U]]), color = 2, label = "")
plot!([-par.ϵ, 0, par.ϵ], collect(solns[3,[:L, :M, :U]]), color = 3, label = "")
scatter!([-par.ϵ, 0, par.ϵ], collect(solns[3,[:L, :M, :U]]), color = 3, label = "")
plot!(zeros(2), [minimum(solns.L), maximum(solns.U)], color = :black, linestyle = :dot, label = "")
plot!([-par.ϵ, par.ϵ], zeros(2), color = :black, linestyle = :dot, label = "")
plot!(size = (300,300))
savefig("figures/analytics/candidate_solns.pdf")


L_htmp = unique(comb_conditions[:,[:L_soln, :M, :U]])
M_htmp = unique(comb_conditions[:,[:L, :U, :M_soln]])
U_htmp = unique(comb_conditions[:,[:L, :M, :U_soln]])

gr()
surface(Ms, Us, Matrix(select!(unstack(L_htmp, :M, :U, :L_soln), Not(:M))))
plot!(xlabel = "M", ylabel = "U", zlabel = "L")
plot!(size = (300,300))
savefig("figures/analytics/L_solns_plane.pdf")

surface(Ls, Us, Matrix(select!(unstack(M_htmp, :L, :U, :M_soln), Not(:L))))
plot!(xlabel = "L", ylabel = "U", zlabel = "M")
plot!(size = (300,300))
savefig("figures/analytics/M_solns_plane.pdf")

surface(Ls, Us, Matrix(select!(unstack(U_htmp, :L, :M, :U_soln), Not(:L))))
plot!(xlabel = "L", ylabel = "M", zlabel = "U")
plot!(size = (300,300))
savefig("figures/analytics/U_solns_plane.pdf")


"""
Define the variables and expectations to keep track of.
All variables which appear as an expectation need to be included here.
Lagged values which appear as state variables need not be included.
"""

@everywhere variables = Symbol.(cat(Vector(options.exogenous),
    Vector(options.endogenous),outputnames(options),
    Vector(options.expectations),Vector(options.auxiliary), dims = 1));
s = DataFrame(zeros(options.N, length(variables)), variables);
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
central[:π] = 0.;
central[:y] = 0.;
central[:Eπ_lead] = 0.;
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
Initialise beliefs by training on (some of the) steady state(s)
"""

@everywhere beliefs = initialise_beliefs(options)
s = initialise_df(s, ss) # gap = 500, steadystate_alt = upper)
s = initialise_df(s, central) # gap = 500, steadystate_alt = upper)
@time beliefs = learn!(beliefs, s, options.N, options, indices, loss)
display(params(beliefs))




"""
Run learning simulation
"""
options.burnin_use_net = true;
options.learning_gap = 10;
options.plotting_gap = 50;
options.window = 499;
options.plot_vars = [:π, :y, :Eπ, :Ey]

par = @set par.ϵ = 10

P = [0.5 0.25 0.25; 0.25 0.5 0.25; 0.25 0.25 0.5]
mc = MarkovChain(P, [-par.ϵ, 0, par.ϵ])
P = [0.4 0.2 0.2 0.2; 0.2 0.4 0.2 0.2; 0.2 0.2 0.4 0.2; 0.2 0.2 0.2 0.4]
mc = MarkovChain(P, [-15, -5, 5, 15])
simulate(mc, 100, init = 2)

gr() # Set GR backend for plots as it's the fastest
# Replace the first shocks with the last ones from the previous time
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_π = simulate(mc, options.N, init = 2)
s.ϵ_y = zeros(options.N)
@time beliefs,s = simulate_learning(options.burnin:options.N, s, beliefs, indices, options)

# Plot simulated time series
pyplot()
nn = 32890;
plot_range = (nn-99):(nn)
plot(layout=(3,1),legend = false,  link = :x)
plot!(s.π[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(s.ϵ_π[plot_range], subplot = 3, ylabel = L"\epsilon_t", yguidefontrotation=-90, xlabel = "Periods")
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

initial_ss = deepcopy(central)
L_est = irf(:ϵ_π, initial_ss, beliefs, shock_period = 2, periods = 2,
	magnitude = -ϵ, persistence = 0., show_plot = false)[2,:]
M_est = irf(:ϵ_π, initial_ss, beliefs, shock_period = 2, periods = 2,
	magnitude = 0., persistence = 0., show_plot = false)[2,:]
U_est = irf(:ϵ_π, initial_ss, beliefs, shock_period = 2, periods = 2,
	magnitude = ϵ, persistence = 0., show_plot = false)[2,:]


"""
Plot phase diagram
"""

pyplot()
plot_ss(plot_points)
initial_ss = deepcopy(central)
starts = [(π=3.0,y=3.0,periods=100,arrows=[10]),
	(π=-3.,y=-3.,periods=100,arrows=[10]),
	(π=0.1,y=0.1,periods=100,arrows=[22, 65]),
	(π=-0.25,y=-0.25,periods=100,arrows=[22, 65]),
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
	plot!(IS_condition.(plot_points),plot_points, label = "IS Curve", color = :green,
		linestyle = :dot)
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

s = initialise_df(s, central);
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
