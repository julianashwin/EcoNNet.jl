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
    expectations = [:Eπ_lead], # [:Eπ,:Ey,:Eπ_lead],
    endogenous = [:π, :y],
    exogenous = [:ϵ_π,:ϵ_y],
    states = [:ϵ_π, :ϵ_y],
	auxiliary = [:r],
    N = 500000, num_nodes = 3, activation = sigmoid, window = 50000,
	burnin = 50000, learning_gap = 50000, plotting_gap = 100)

@everywhere beliefs = initialise_beliefs(options)
params(beliefs)

"""
Define the parameters as a Named Tuple.
"""
## Model parameters
@everywhere par = (β = 1, κ = 1, η = 0, σ = 1, ζ = 0.2,
	ϕ_π = 0, π_star = 1.0, α = 5.0, ϵ = [-2.0, -1.0, -1/3, 1/3, 1.0, 2.0],
	p = 0.75);
# Shock transition matrix
P = [par.p (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 (1-par.p)/5;
    (1-par.p)/5 par.p (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 (1-par.p)/5;
    (1-par.p)/5 (1-par.p)/5 par.p (1-par.p)/5 (1-par.p)/5 (1-par.p)/5;
    (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 par.p (1-par.p)/5 (1-par.p)/5;
    (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 par.p (1-par.p)/5;
    (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 (1-par.p)/5 par.p]
heatmap(P)
# Define as a MarkovChain, using the routines from QuantEcon.jl
ϵ_π_mc = MarkovChain(P, par.ϵ)
simulate(ϵ_π_mc, 100, init = 3)
    

"""
State the equilibrium conditions of the model as a function which returns
    a vector of zeros
"""

@everywhere function equilibrium_conditions_fast(F::Array{Float64,1},
    x::Array{Float64,1},states_input::Array{Float64,1},predictions_input::Array{Float64,1})
    # Manually unpack the states
    y_lag::Float64 = 0.#states_input[1]
    ϵ_π::Float64 = states_input[1]
	ϵ_y::Float64 = states_input[2]
    # and the predictions
    Eπ_lead::Float64 = predictions_input[1]
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


# Are there multiple 2,2,2 MSV solutions?
σ̃ = par.σ*par.κ;
par = @set par.p = 0.75
function C_cond(par, D) 
    σ̃ = par.σ*par.κ;
    C = (1/(1 - (1+σ̃)*(par.p)))*((1 + σ̃)*((1 - par.p)/5)*D + par.ϵ[3])
end
function D_cond(par, C) 
    σ̃ = par.σ*par.κ;
    D = (1/(1 - (1+σ̃)*(par.p)))*((1 + σ̃)*((1 - par.p)/5)*C + par.ϵ[4])
end
plot(C_cond.([par], -1.0:0.01:1.0), -1.0:0.01:1.0, xlabel = "C", ylabel = "D",
    label = L"C = \frac{1}{1 - (1 + \tilde{\sigma})p}\left[  \frac{(1 + \tilde{\sigma})(1-p)}{5} D  - \frac{1}{3}  \right]")
plot!(-1.0:0.01:1.0, D_cond.([par], -1.0:0.01:1.0), 
    label = L"D = \frac{1}{1 - (1 + \tilde{\sigma})p}\left[  \frac{(1 + \tilde{\sigma})(1-p)}{5} C  + \frac{1}{3}  \right]")
plot!(title = "p = "*string(par.p))




"""
Identify potential REEs numerically
"""
Z = zeros(6)
function check_REE(par, Z)
    A = Z[1]; B = Z[2]; C = Z[3]; D = Z[4]; E = Z[5]; F = Z[6]; 
    σ̃ = par.σ*par.κ;
    z_ = zeros(6)
    # A condition 
    if A < -par.π_star
        z_[1] = (1 + σ̃)*((par.p)*A + ((1 - par.p)/5)*(B + C + D + E + F )) - par.α*σ̃*(A + par.π_star) + par.ϵ[1]
    elseif (A >= -par.π_star) & (A <= par.π_star) 
        z_[1] = (1 + σ̃)*((par.p)*A + ((1 - par.p)/5)*(B + C + D + E + F )) + par.ϵ[1]
    elseif A > par.π_star
        z_[1] = (1 + σ̃)*((par.p)*A + ((1 - par.p)/5)*(B + C + D + E + F )) - par.α*σ̃*(A - par.π_star) + par.ϵ[1]
    end 
    # B condition 
    if B < -par.π_star
        z_[2] = (1 + σ̃)*((par.p)*B + ((1 - par.p)/5)*(A + C + D + E + F )) - par.α*σ̃*(B + par.π_star) + par.ϵ[2]
    elseif (B >= -par.π_star) & (B <= par.π_star) 
        z_[2] = (1 + σ̃)*((par.p)*B + ((1 - par.p)/5)*(A + C + D + E + F )) + par.ϵ[2]
    elseif B > par.π_star
        z_[2] = (1 + σ̃)*((par.p)*B + ((1 - par.p)/5)*(A + C + D + E + F )) - par.α*σ̃*(B - par.π_star) + par.ϵ[2]
    end 
    # C condition
    if C < -par.π_star
        z_[3] = (1 + σ̃)*((par.p)*C + ((1 - par.p)/5)*(A + B + D + E + F )) - par.α*σ̃*(C + par.π_star) + par.ϵ[3]
    elseif (C >= -par.π_star) & (C <= par.π_star) 
        z_[3] = (1 + σ̃)*((par.p)*C + ((1 - par.p)/5)*(A + B + D + E + F )) + par.ϵ[3]
    elseif C > par.π_star
        z_[3] = (1 + σ̃)*((par.p)*C + ((1 - par.p)/5)*(A + B + D + E + F )) - par.α*σ̃*(C - par.π_star) + par.ϵ[3]
    end 
    # D condition 
    if D < -par.π_star
        z_[4] = (1 + σ̃)*((par.p)*D + ((1 - par.p)/5)*(A + B + C + E + F )) - par.α*σ̃*(D + par.π_star) + par.ϵ[4]
    elseif (D >= -par.π_star) & (D <= par.π_star) 
        z_[4] = (1 + σ̃)*((par.p)*D + ((1 - par.p)/5)*(A + B + C + E + F )) + par.ϵ[4]
    elseif D > par.π_star
        z_[4] = (1 + σ̃)*((par.p)*D + ((1 - par.p)/5)*(A + B + C + E + F )) - par.α*σ̃*(D - par.π_star) + par.ϵ[4]
    end 
    # E condition
    if E < -par.π_star
        z_[5] = (1 + σ̃)*((par.p)*E + ((1 - par.p)/5)*(A + B + C + D + F )) - par.α*σ̃*(E + par.π_star) + par.ϵ[5]
    elseif (E >= -par.π_star) & (E <= par.π_star) 
        z_[5] = (1 + σ̃)*((par.p)*E + ((1 - par.p)/5)*(A + B + C + D + F )) + par.ϵ[5]
    elseif E > par.π_star
        z_[5] = (1 + σ̃)*((par.p)*E + ((1 - par.p)/5)*(A + B + C + D + F )) - par.α*σ̃*(E - par.π_star) + par.ϵ[5]
    end 
    # F condition
    if F < -par.π_star
        z_[6] = (1 + σ̃)*((par.p)*F + ((1 - par.p)/5)*(A + B + C + D + E )) - par.α*σ̃*(F + par.π_star) + par.ϵ[6]
    elseif (F >= -par.π_star) & (F <= par.π_star) 
        z_[6] = (1 + σ̃)*((par.p)*F + ((1 - par.p)/5)*(A + B + C + D + E )) + par.ϵ[6]
    elseif F > par.π_star
        z_[6] = (1 + σ̃)*((par.p)*F + ((1 - par.p)/5)*(A + B + C + D + E )) - par.α*σ̃*(F - par.π_star) + par.ϵ[6]
    end 
    
    return z_
end

# Numerically solve for different REEs 
function REE_solve(z; par = par)
    z_ = check_REE(par, z)
    out = z_ - z
    return out
end

par = @set par.p = 0.75
par = @set par.π_star = 1.0
soln0 = nlsolve(REE_solve, zeros(length(par.ϵ)))
#soln1 = nlsolve(REE_solve, π_msv)
soln1 = nlsolve(REE_solve, par.π_star.*[-1.3, -1.3, -0.3, 0.3, 1.3, 1.3])
soln2 = nlsolve(REE_solve, par.π_star.*[-1.3, -1.3, -1.3, 1.3, 1.3, 1.3])
soln3 = nlsolve(REE_solve, par.π_star.*[-1.3, -1.3, 1.3, 1.3, 1.3, 1.3])
soln4 = nlsolve(REE_solve, par.π_star.*[-1.3, -1.3, -1.3, -1.3, 1.3, 1.3])
soln5 = nlsolve(REE_solve, par.π_star.*[-1.3, 1.3, 1.3, 1.3, 1.3, 1.3])
soln6 = nlsolve(REE_solve, par.π_star.*[-1.3, -1.3, -1.3, -1.3, -1.3, 1.3])
soln7 = nlsolve(REE_solve, par.π_star.*[1.3, 1.3, 1.3, 1.3, 1.3, 1.3])
soln8 = nlsolve(REE_solve, 2.0.*par.π_star.*[-1.3, -1.3, -1.3, -1.3, -1.3, -1.3])


pyplot()
plot(layout = (1), legend = :outerright, link = :y)
plot!(xlabel = L"\epsilon_t", ylabel = L"\pi_t", 
    title = L"p = "*string(par.p)*L", \pi^* = "*string(par.π_star))
plot!(par.ϵ, par.π_star*ones(len((par.ϵ))), fillrange=[-par.π_star*ones(len(par.ϵ))], fillalpha = 0.5,
		 color = :paleturquoise1, label = L"|\pi_t| < \pi^*")
if sum(abs.(soln0.zero - check_REE(par, soln0.zero))) < 1e-6
    scatter!(par.ϵ, soln0.zero, label = "REE 0 "*string(round.(soln0.zero, digits = 2)) , color = 0)
    plot!(par.ϵ, soln0.zero, label = false, color = 0)
end
if sum(abs.(soln1.zero - check_REE(par, soln1.zero))) < 1e-6
    scatter!(par.ϵ, soln1.zero, label = "REE 1 "*string(round.(soln1.zero, digits = 2)) , color = 1)
    plot!(par.ϵ, soln1.zero, label = false, color = 1)
end
if sum(abs.(soln2.zero - check_REE(par, soln2.zero))) < 1e-6
    scatter!(par.ϵ, soln2.zero, label = "REE 2 "*string(round.(soln2.zero, digits = 2)), color = 2)
    plot!(par.ϵ, soln2.zero, label = false, color = 2)
end
if sum(abs.(soln3.zero - check_REE(par, soln3.zero))) < 1e-6
    scatter!(par.ϵ, soln3.zero, label = "REE 3 "*string(round.(soln3.zero, digits = 2)), color = 3)
    plot!(par.ϵ, soln3.zero, label = false, color = 3)
end
if sum(abs.(soln4.zero - check_REE(par, soln4.zero))) < 1e-6
    scatter!(par.ϵ, soln4.zero, label = "REE 4 "*string(round.(soln4.zero, digits = 2)), color = 4)
    plot!(par.ϵ, soln4.zero, label = false, color = 4)
end
if sum(abs.(soln5.zero - check_REE(par, soln5.zero))) < 1e-6
    scatter!(par.ϵ, soln5.zero, label = "REE 5 "*string(round.(soln5.zero, digits = 2)), color = 5)
    plot!(par.ϵ, soln5.zero, label = false, color = 5)
end
if sum(abs.(soln6.zero - check_REE(par, soln6.zero))) < 1e-6
    scatter!(par.ϵ, soln6.zero, label = "REE 6 "*string(round.(soln6.zero, digits = 2)), color = 6)
    plot!(par.ϵ, soln6.zero, label = false, color = 6)
end
if sum(abs.(soln7.zero - check_REE(par, soln7.zero))) < 1e-6
    scatter!(par.ϵ, soln7.zero, label = "REE 7 "*string(round.(soln7.zero, digits = 2)), color = 7)
    plot!(par.ϵ, soln7.zero, label = false, color = 7)
end
if sum(abs.(soln8.zero - check_REE(par, soln8.zero))) < 1e-6
    scatter!(par.ϵ, soln8.zero, label = "REE 8 "*string(round.(soln8.zero, digits = 2)), color = 8)
    plot!(par.ϵ, soln8.zero, label = false, color = 8)
end
soln_plt = plot!(size = (1000,400))
plot(soln_plt)
savefig("figures/analytics/candidate_REEs.pdf")




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
Initialise beliefs by training on  (potentially peturbed) REE solution
"""
# Potential as and bs for (2,2,2) solution
a1 = par.α*par.σ*par.κ*par.π_star
a2 = (6*(1 + par.σ*par.κ)*par.p - 5*par.σ*par.κ)*par.σ*par.κ/5
a3 = -(6*(1 + par.σ*par.κ)*par.p - 5*par.σ*par.κ)*par.σ*par.κ/5
b11 = 1.0
b12 = 1.0
b13 = 1.0
b21 = - (6*par.p - 5) / (6*(1 + par.σ*par.κ)*par.p - 5*par.σ*par.κ*(par.α - 1))
b22 = 5*par.α*(6*par.p - 5) / ((6*(1 + par.σ*par.κ)*par.p - 5*par.σ*par.κ) * 
    (6*(1 + par.σ*par.κ)*par.p - 5*par.σ*par.κ*(par.α - 1)))
b23 = - 5*par.α*(6*par.p - 5) / ((6*(1 + par.σ*par.κ)*par.p - 5*par.σ*par.κ) * 
    (6*(1 + par.σ*par.κ)*par.p - 5*par.σ*par.κ*(par.α - 1)))

if (false)
a1 = par.α*par.σ*par.κ*par.π_star
a2 = (6*(1 + par.σ*par.κ)*(1-par.p) - 5*par.σ*par.κ)*par.σ*par.κ/5
a3 = -(6*(1 + par.σ*par.κ)*(1-par.p) - 5*par.σ*par.κ)*par.σ*par.κ/5
b21 = - (6*(1-par.p) - 5) / (6*(1 + par.σ*par.κ)*(1-par.p) - 5*par.σ*par.κ*(par.α - 1))
b22 = 5*par.α*(6*(1-par.p) - 5) / ((6*(1 + par.σ*par.κ)*(1-par.p) - 5*par.σ*par.κ) * 
    (6*(1 + par.σ*par.κ)*(1-par.p) - 5*par.σ*par.κ*(par.α - 1)))
b23 = - 5*par.α*(6*(1-par.p) - 5) / ((6*(1 + par.σ*par.κ)*(1-par.p) - 5*par.σ*par.κ) * 
    (6*(1 + par.σ*par.κ)*(1-par.p) - 5*par.σ*par.κ*(par.α - 1)))
end


params(beliefs)[1] .= hcat(Float32.([1.0, 1.0, 1.0])) # b1
params(beliefs)[2] .= Float32.([a1, a2, a3]) # a1
params(beliefs)[3] .= Float32.([b21 b22 b23]) # b2
plot(par.ϵ, vcat(beliefs([par.ϵ[1]]), beliefs([par.ϵ[2]]), beliefs([par.ϵ[3]]), beliefs([par.ϵ[4]]), 
    beliefs([par.ϵ[5]]), beliefs([par.ϵ[6]])))



b21.*(a1 .+ b11.*par.ϵ)
max.(a2 .+ b12.*par.ϵ, [0])
max.(a3 .+ b13.*par.ϵ, [0])
inpt = par.ϵ[1]
b21*(a1 + b11*inpt) + b22*max(a2 + b12*inpt, 0) + b23*max(a3 + b13*inpt, 0)

beliefs([par.ϵ[1]])

function init_eqlm(par, s, soln)
    for ii in 1:length(par.ϵ)
        exp_inf = par.p*soln[ii] + (1-par.p)/5*sum(soln[Not(ii)])
        s.π[(s.ϵ_π .== par.ϵ[ii])] .= soln[ii]
        s.Eπ_lead[(s.ϵ_π .== par.ϵ[ii])] .= exp_inf
    end
    s.r .= 0.0
    s.r[s.π .> par.π_star] = par.ϕ_π.*s.π[s.π .> par.π_star] .+ par.α.*(s.π[s.π .> par.π_star] .- par.π_star)
    s.r[s.π .< -par.π_star] = par.ϕ_π.*s.π[s.π .< -par.π_star] .+ par.α.*(s.π[s.π .< -par.π_star] .+ par.π_star)
	s.y = par.σ.*(s.r .- s.Eπ_lead) .+ s.ϵ_y
    return s
end
function learned_eqlm(par, beliefs)
    π_eqlm = zeros(length(par.ϵ))
    for ii in 1:length(par.ϵ)
        inputs = [par.ϵ[ii]]::Array{Float64,1}
        predictions = predict!(inputs, beliefs)
        states = [par.ϵ[ii], 0.0]::Array{Float64,1}
        outcomes = step_fast!(cat([0.; 0.], states, predictions, dims = 1), options)
        π_eqlm[ii] = outcomes[1]
    end
    return π_eqlm
end
# Only need to initialise ϵ_π and π
options.window = 50000
options.num_nodes = 3
options.activation = relu
options.optim = Descent(0.4)
#options.optim = ADAM()
function loss(x, y)
	Flux.mae(beliefs(x), y)
end
function loss(x, y)
	mse(beliefs(x), y)
end


@everywhere beliefs = initialise_beliefs(options)
params(beliefs)[1] .= hcat(Float32.([1.0, 1.0, 1.0])) # b1
params(beliefs)[2] .= Float32.([2., 0., 0.]) # a1
params(beliefs)[3] .= Float32.([1.0 1.0 1.0]) # b2

params(beliefs)[1] .= hcat(Float32.([1.0, 1.0, 1.0])) # b1
params(beliefs)[2] .= Float32.([a1, a2, a3]) # a1
params(beliefs)[3] .= Float32.([b21 b22 b23]) # b2
plot(par.ϵ, vcat(beliefs([par.ϵ[1]]), beliefs([par.ϵ[2]]), beliefs([par.ϵ[3]]), beliefs([par.ϵ[4]]), 
    beliefs([par.ϵ[5]]), beliefs([par.ϵ[6]])))

π_learned = learned_eqlm(par, beliefs)   

init_soln = soln2.zero
s.ϵ_π = simulate(ϵ_π_mc, options.N)
s.ϵ_y = randn(options.N)
s = init_eqlm(par, s, init_soln)
scatter(s.ϵ_π[1:500], s.π[1:500], label = "REE", markershape = :xcross, markersize = 6, legend = :outerright)
scatter!(s.ϵ_π[1:500], s.Eπ_lead[1:500], label = "Exp", legend = :outerright)

# Check for REE with DHM stat
pvals = DHM_test(s, 5000:500:500000, 500, hvars = [:ϵ_π], include_const = true);



# Train beliefs on this initialisation
@time beliefs = learn!(beliefs, s, options.N-1, options, indices, loss)
π_learned = learned_eqlm(par, beliefs)   
scatter(s.ϵ_π[1:500], s.π[1:500], label = "REE", markershape = :xcross, markersize = 6, legend = :outerright)
scatter!(par.ϵ, π_learned, label = "Beliefs")






"""
Run learning simulation
"""
options.burnin_use_net = true;
options.learning_gap = 50000;
options.plotting_gap = 50000;
options.window = 50000;
options.plot_vars = [:π, :y, :Eπ_lead]


# Replace the first shocks with the last ones from the previous time
s[1:options.burnin,:] = s[(options.N-options.burnin+1):options.N,:]
s.ϵ_π[(options.burnin+1):options.N] = simulate(ϵ_π_mc, options.N-options.burnin, init = 2)
s.ϵ_y = zeros(options.N)
gr() # Set GR backend for plots as it's the fastest
@time beliefs,s = simulate_learning(options.burnin+1:options.N, s, beliefs, indices, options)

π_learned = learned_eqlm(par, beliefs)   
display(soln_plt)
plot(soln_plt)
scatter!(par.ϵ, π_learned, label = "Learnable", markershape = :xcross, markersize = 6, 
    color = :black, legend = :outerright)

pvals = DHM_test(s, 5000:500:500000, 200, hvars = [:ϵ_π], include_const = true);


# Plot the equilibrium
pyplot()
plot_df = unique(s[(nn-1999):(nn),:])
plot(layout = (2,1), legend = :topleft, link = :y)
scatter!(plot_df.ϵ_π, plot_df.π, markerstrokewidth = 0, label = "Learnable")
plot!(xlabel = L"\epsilon_t", ylabel = L"\pi_t")
plot!(par.ϵ, par.π_star*ones(len((par.ϵ))), fillrange=[-par.π_star*ones(len(par.ϵ))], fillalpha = 0.5,
		 color = :paleturquoise1, label = L"\pi_t < \pi^*")
scatter!(par.ϵ, π_msv, markershape = :xcross, color = :green, markerstrokewidth = 0, label = "MSV solution")
scatter!(plot_df.ϵ_π, plot_df.Eπ_lead, markerstrokewidth = 0, label = "Learnable", subplot = 2)
plot!(xlabel = L"\epsilon_t", ylabel = L"E \pi_{t+1}", legend = :none, subplot = 2)
plot!(par.ϵ, par.π_star*ones(len((par.ϵ))), fillrange=[-par.π_star*ones(len(par.ϵ))], fillalpha = 0.5,
		 color = :paleturquoise1, label = L"\pi_t < \pi^*", subplot = 2)
scatter!(par.ϵ, π_lead_msv, markershape = :xcross, color = :green, markerstrokewidth = 0, label = "MSV solution", subplot = 2)
plot!(size = (600,450))


"""
Check for REE
    inputs are [ϵ_π]
    predictions are [Eπ_lead]
    states are [ϵ_π, ϵ_y]
    outcomes are [π_t, y_t]

"""



check_REE(par, π_eqlm)




# Conduct DHM test
range = 40000:50000 .-1000
errors = s.π[range.-1] - s.Eπ_lead[range]
mean(errors)


# Plot simulated time series
pyplot()
nn = 50000;
plot_range = (nn-99):(nn)
plot(layout=(3,1),legend = false,  link = :x)
plot!(s.π[plot_range], subplot = 1, ylabel = L"\pi_t", yguidefontrotation=-90)
plot!(s.y[plot_range], subplot = 2, ylabel = L"y_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(s.ϵ_π[plot_range], subplot = 3, ylabel = L"\epsilon_t", yguidefontrotation=-90, xlabel = "Periods")
plot!(size = (600,300))


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
Compare PLM and ALM
"""
s.Eπ_lead_error = vcat((s.Eπ_lead[1:(options.N-1)]-s.π[2:(options.N)]),[0.])
plot(s.Eπ_lead_error)

nn=49000
plot_range = (nn-9999):(nn)
scatter(s.ϵ_π[plot_range], s.Eπ_lead_error)


avg_error = zeros(len(par.ϵ))
for ii in 1:len(par.ϵ)
	short_df = s[plot_range,:]
	instances = Array(short_df[:,:ϵ_π].==par.ϵ[ii])
	errors = short_df[instances,:Eπ_lead_error]
	avg_error[ii] = mean(errors)
	if isnan(avg_error[ii])
		avg_error[ii] = 0.0
	end
end
plot(par.ϵ, avg_error, xlabel= L"\epsilon_{\pi,t}", ylabel= L"Avg. [E \pi_{t+1} - \pi_{t+1}]",
	legend=:false,yguidefontrotation=0)
plot!(size=(600,400))
savefig("figures/heatmap_errors_illus_sim.pdf")



Rsq_π = 1 - var(heatmap_df.π-heatmap_df.Eπ)/var(heatmap_df.π)
Rsq_y = 1 - var(heatmap_df.y-heatmap_df.Ey)/var(heatmap_df.y)







"""
Characterise the REE
"""
## MSV solution with two in middle (from Eq 21 in Martin's note)
σ̃ = par.σ*par.κ
π_msv = zeros(6)
π_msv[1] = (5*(par.ϵ[1] - par.α*σ̃))/(6*(1 + σ̃)*par.p + 5*σ̃*(par.α - 1))
π_msv[2] = (5*(par.ϵ[2] - par.α*σ̃))/(6*(1 + σ̃)*par.p + 5*σ̃*(par.α - 1))
π_msv[3] = (5*par.ϵ[3])/(6*(1 + σ̃)*par.p - 5*σ̃)
π_msv[4] = (5*par.ϵ[4])/(6*(1 + σ̃)*par.p - 5*σ̃)
π_msv[5] = (5*(par.ϵ[5] + par.α*σ̃))/(6*(1 + σ̃)*par.p + 5*σ̃*(par.α - 1))
π_msv[6] = (5*(par.ϵ[6] + par.α*σ̃))/(6*(1 + σ̃)*par.p + 5*σ̃*(par.α - 1))

π_lead_msv = zeros(6)
π_lead_msv[1] = -((6*par.p - 5)*(par.ϵ[1] - par.α*σ̃))/(6*(1 + σ̃)*par.p + 5*σ̃*(par.α - 1))
π_lead_msv[2] = -((6*par.p - 5)*(par.ϵ[2] - par.α*σ̃))/(6*(1 + σ̃)*par.p + 5*σ̃*(par.α - 1))
π_lead_msv[3] = -((6*par.p - 5)*par.ϵ[3])/(6*(1 + σ̃)*par.p - 5*σ̃)
π_lead_msv[4] = -((6*par.p - 5)*par.ϵ[4])/(6*(1 + σ̃)*par.p - 5*σ̃)
π_lead_msv[5] = -((6*par.p - 5)*(par.ϵ[5] + par.α*σ̃))/(6*(1 + σ̃)*par.p + 5*σ̃*(par.α - 1))
π_lead_msv[6] = -((6*par.p - 5)*(par.ϵ[6] + par.α*σ̃))/(6*(1 + σ̃)*par.p + 5*σ̃*(par.α - 1))

π_lead_msv[1] 
sum(P[1,:].*π_msv)

scatter(π_msv)
scatter!(check_REE(par, π_msv))

# Plot the first MSV REE
pyplot()
plot_df = unique(s[(nn-1999):(nn),:])
plot(layout = (2,1), legend = :bottomright, link = :y)
plot!(xlabel = L"\epsilon_t", ylabel = L"\pi_t")
plot!(par.ϵ, par.π_star*ones(len((par.ϵ))), fillrange=[-par.π_star*ones(len(par.ϵ))], fillalpha = 0.5,
            color = :paleturquoise1, label = L"\pi_t < \pi^*")
scatter!(par.ϵ, π_msv, markershape = :xcross, color = :green, markerstrokewidth = 0, label = "MSV solution")
plot!(par.ϵ, par.π_star*ones(len((par.ϵ))), fillrange=[-par.π_star*ones(len(par.ϵ))], fillalpha = 0.5,
            color = :paleturquoise1, label = L"\pi_t < \pi^*", subplot = 2)
plot!(xlabel = L"\epsilon_t", ylabel = L"E \pi_{t+1}", legend = :none, subplot = 2)
scatter!(par.ϵ, π_lead_msv, markershape = :xcross, color = :green, markerstrokewidth = 0, label = "MSV solution", subplot = 2)
plot!(size = (600,450))