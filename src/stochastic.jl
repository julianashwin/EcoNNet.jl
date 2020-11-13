"""
Functions to deal with the exogenous stochastic shocks
"""


"""
A function that generates n draws from an AR process with persistence
	ρ and noise parameterised by σ drawn from dist distribution.
"""
function simulate_ar(ρ::Float64, σ::Float64, n::Int64, dist)
	noise::Array{Float64,1} = zeros(n)
	if len(dist) == n
		noise = dist
	else
		noise = σ*rand(dist, n)
	end
	y::Array{Float64,1} = [0.0 for ii = 1:n]
	y[1] = noise[1]
	for ii in 2:(n)
		y[ii] = ρ*y[ii-1] + noise[ii]
	end
	return y
end


"""
A function that creates rare large shocks, of alternating sign.
    gap is the gap between shocks
    mag is the magnitude of the shock (in absolute terms)
    start is the first period in which the shock hits
    N is the total length of the shock series
"""
function alternating_shocks(N::Int64; gap::Int64=100, mag::Float64 = 1.0,
		start::Int64=1)
    noise::Array{Float64,1} = zeros(N)
    noise[start:(2*gap):N] = mag*ones(len(start:(2*gap):N))
    noise[(start+gap):(2*gap):N] = -mag*ones(len((start+gap):(2*gap):N))
    return noise
end


"""
A function that returns support and transition kernel for an N-state
	Markov chain Rouwenhorst approximation of an AR(1) process
	p:
	q:
	m: must be a positive scalar
	N: number of states, must be a positive integer
"""
function Rouwenhorst(p::Float64,q::Float64,ψ::Float64,N::Int)

	vtest = [p, q, ψ, N]
	@assert (size(vtest,1) == 4) "Rouwenhorst(): argument dimensions incorrect -- must all be scalars."
	@assert (p >= 0 && p < 1) "Rouwenhorst(): p must be between zero and 1."
	@assert (q >= 0 && q < 1) "Rouwenhorst(): q must be between zero and 1."
	@assert (ψ > 0) "Rouwenhorst(): m must be positive."
 	@assert (N > 2) "Rouwenhorst(): N must be an integer greater than 1."

	global mTheta = [p (1-p); 1-q q]
	for ii = 2:1:(N-1)
		temp = p*[mTheta zeros(ii,1); zeros(1, ii+1)] +
			(1-p) * [zeros(ii,1) mTheta; zeros(1, ii+1)] +
		  	(1-q) * [zeros(1, ii+1); mTheta zeros(ii,1)] +
			q * [zeros(1, ii+1); zeros(ii,1) mTheta]
		mTheta=temp
		mTheta[2:ii,:] = mTheta[2:ii,:]/2;
	end
	vm = -ψ : 2*ψ/(N-1) : ψ
	vm = Array(vm)
	if len(vm) < N
		vm = cat(vm, ψ, dims = 1)
	end

	return vm, mTheta
end




"""
Function that creates a transition kernel from a dataframe with transition probabilities
"""
function create_kernel(df::DataFrame, poss_states::Int64)

	# Create the empty transition probability matrix
	CP::Array{Float64,2} = zeros(poss_states, poss_states)

	# Pick out the relevant variables from df
	trans_probs::Array{Real,2} = Array{Real,2}(df_grid_new[:,[:state_t,:state_tp1,
		:trans_prob]])

	# Populate the CP matrix with the transition probabilities
	prog = Progress(size(trans_probs,1), dt = 1, desc="Populating CP matrix: ")
	for st in 1:size(trans_probs,1)
		CP[trans_probs[st,1],trans_probs[st,2]] += trans_probs[st,3]
		next!(prog)
	end

	return CP
end



"""
Compute stationary distribution from stochastic matrix
"""
function MC_stationary(A::Array{Float64,2}, max_iter::Int64 = 100)
    n_states::Int64 = size(A,1)
    @assert round(sum(A),digits = 10) == n_states "Must be valid stochastic matrix (rows sum to one)"

    prog = Progress(6, dt = 1, desc="Computing marginal probability eigenvector: ")
    next!(prog)
    # Identify which index to set to one
    idx::Int64 = 0
    for ii in 1:max_iter
        idx = sample(1:n_states)
        if sum(A[:,idx]) > (1/n_states)
            break
        end
    end
    next!(prog)

    idxs::Array{Int64,1} = vcat(Array(1:(idx-1)),Array((idx+1):n_states))
    numer = reshape(-A[idx,idxs],(n_states-1,1))
    next!(prog)
    A2::Array{Float64,2} = Matrix(transpose(A[idxs,idxs]))
    next!(prog)
    denom::Array{Float64,2} = A2-I(n_states-1)
    next!(prog)
    x::Array{Float64,1}=ones(n_states)
    x[idxs] = \(denom,numer)
    next!(prog)
    x = x./sum(x)
    next!(prog)
    return x
end

"""
Faster version of MC_stationary that uses the fact that some columns are all zeros
"""
function MC_stationary_fast(A::Array{Float64,2}, cutoff::Float64 = 0.)
    n_states::Int64 = size(A,1)
    @assert round(sum(A),digits = 10) == n_states "Must be valid stochastic matrix (rows sum to one)"

    keep_cols::Array{Int64,1} = (1:n_states)[vec(sum(A,dims=1)).>cutoff]
    A_short::Array{Float64,2} = A[keep_cols,keep_cols]
    x_short::Array{Float64,1} = MC_stationary(A_short)
    x::Array{Float64,1} = zeros(n_states)
    x[keep_cols] = x_short

    return x
end
