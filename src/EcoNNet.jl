"""
General convention on ordering variables:
	1) Order in chronological order (i.e. lagged, current, lead)
	2) Order alphabetically within each period
"""

# Import dependencies
using Distributions, Parameters, LinearAlgebra, QuantEcon
using Flux, Distributed, SharedArrays, JLD, NLsolve, ProgressMeter
using Flux: throttle, params, mse, glorot_uniform, @epochs
using IterTools: ncycle, NCycle
using NLsolve: SolverResults
using DataFrames, Random, Atom, Juno, DataFramesMeta
using Plots, Plots.PlotMeasures



"""
Mutable structure for options
"""
@with_kw mutable struct EcoNNetOptions
    N::Int64 = 100000;
    window::Int64 = 5000;
    NNet::Bool = true;
    init_weights = glorot_uniform;
    optim = ADAM;
    activation = σ;
    hidden_layers::Int64 = 1;
    train_split::Array{Float64,1} = [1.0, 0.0, 0.0];
    num_nodes::Int64 = 24;
    max_iter::Int64 = 10;
    infoset::Array{Symbol,1};
    expectations::Array{Symbol,1};
    outputs::Array{Symbol,1} = Symbol.([]);
    endogenous::Array{Symbol,1};
    exogenous::Array{Symbol,1};
    states::Array{Symbol,1};
    auxiliary::Array{Symbol,1} = Symbol.([]);
	burnin::Int64 = 6000;
	burnin_use_net::Bool = false;
	learning_gap::Int64 = 100;
	plotting_gap::Int64 = 100;
	plotting_window::Int64 = 999;
	show_plots::Bool = true;
	plot_vars::Array{Symbol,1} = Symbol.([]);
end

"""
Mutable struct with indices of the different types of variable
"""
@with_kw mutable struct EcoNNetIndices
	constant::Int64 = 0;
	infoindex_lagged::Array{Int64,1};
	infonames_lagged::Array{Symbol,1};
	infoindex_current::Array{Int64,1};
	infonames_current::Array{Symbol,1};
	outputindex_current::Array{Int64,1};
	outputnames_current::Array{Symbol,1};
	outputindex_lead::Array{Int64,1};
	outputnames_lead::Array{Symbol,1};
	expectindex_current::Array{Int64,1};
	expectnames_current::Array{Symbol,1};
	expectindex_lead::Array{Int64,1};
	expectnames_lead::Array{Symbol,1};
	expectindex_all::Array{Int64,1};
	expectnames_all::Array{Symbol,1};
	stateindex_lagged::Array{Int64,1};
	statenames_lagged::Array{Symbol,1};
	stateindex_current::Array{Int64,1} ;
	statenames_current::Array{Symbol,1};
	stateindex_all::Array{Int64,1};
	statenames_all::Array{Symbol,1};
	endogindex::Array{Int64,1};
	endognames::Array{Symbol,1};
	auxindex::Array{Int64,1};
	auxnames::Array{Symbol,1};
end

"""
A shortened name for the length function to make code more compact
"""
function len(x)
	length(x)
end

"""
A function that identifies the position of elements in an array while preserving their order
"""
function find_index(element_names::Array{Symbol,1}, variables::Array{Symbol,1})
	element_index::Array{Int64,1} = fill(0, len(element_names))
	for ii in 1:len(element_names)
		element_index[ii] = findall(x-> x in [element_names[ii]], variables)[1];
	end
	return element_index
end

"""
A function that creates output names from the stated expectations,
	if they are different from just the endogenous variables
"""
function outputnames(options)
	expectations = options.expectations
	outputs = Vector{Symbol}(undef,len(expectations))
	for ii in 1:len(expectations)
		temp = split(string(expectations[ii]), "_")[1]
		temp = split(temp, "E")[2]
		outputs[ii] = Symbol(temp)
	end
	filter!(x->!(x in options.endogenous),outputs)
	filter!(x->!(x in options.exogenous),outputs)
	return outputs
end


"""
A function that indentifies the indices in s for the information set and expectations
"""
function extract_indices(options::EcoNNetOptions, variables::Array{Symbol,1})
	infonames_lagged = Array{Symbol,1}(undef,0);
	infonames_current = Array{Symbol,1}(undef,0);
	infonames_all = Array{Symbol,1}(undef,0);
	outputnames_current = Array{Symbol,1}(undef,0);
	outputnames_lead = Array{Symbol,1}(undef,0);
	outputnames_all = Array{Symbol,1}(undef,0);
	statenames_lagged = Array{Symbol,1}(undef,0);
	statenames_current = Array{Symbol,1}(undef,0);
	statenames_all = Array{Symbol,1}(undef,0);
	expectnames_current = Array{Symbol,1}(undef,0);
	expectnames_lead = Array{Symbol,1}(undef,0);
	expectnames_all = Array{Symbol,1}(undef,0);

	constant = 0;

	auxnames::Array{Symbol,1} = options.auxiliary;
	endognames::Array{Symbol,1} = options.endogenous;

	# indices for information set, split into lagged and current values
    for ii in 1:len(options.infoset)
        if occursin("lag", string(options.infoset[ii]))
            push!(infonames_lagged, Symbol(split(string(options.infoset[ii]) ,"_lag")[1]))
        elseif occursin("constant", string(options.infoset[ii]))
			constant = 1
		else
            push!(infonames_current, options.infoset[ii]);
        end
    end
	infonames_all = cat(infonames_lagged,infonames_current, dims = 1)
	# indices for outputs, split into next and this period
	for ii in 1:len(options.expectations)
        if occursin("lead", string(options.expectations[ii]))
			temp = split(string(options.expectations[ii]) ,"_lead")[1]
			temp = split(temp ,"E")[2]
            push!(outputnames_lead, Symbol(temp))
        else
            push!(outputnames_current, Symbol(split(string(options.expectations[ii]) ,"E")[2]));
        end
    end
	outputnames_all = cat(outputnames_current, outputnames_lead, dims = 1)
	# indices for expectations, split into next and this period
	for ii in 1:len(options.expectations)
        if occursin("lead", string(options.expectations[ii]))
			push!(expectnames_lead, options.expectations[ii])
        else
            push!(expectnames_current, options.expectations[ii]);
        end
    end
	# The expectations are ordered current then lead, for populating the s DataFrame
	expectnames_all = cat(expectnames_current, expectnames_lead, dims = 1)

	# indices for state variables
	for ii in 1:len(options.states)
        if occursin("lag", string(options.states[ii])) && !occursin("lag2", string(options.states[ii]))
            push!(statenames_lagged, Symbol(split(string(options.states[ii]) ,"_lag")[1]))
        else
            push!(statenames_current, options.states[ii]);
        end
    end
	statenames_all = cat(statenames_lagged, statenames_current, dims = 1)

	infoindex_lagged::Array{Int64,1} = find_index(infonames_lagged, variables);
	infoindex_current::Array{Int64,1} = find_index(infonames_current, variables);
	infoindex_all::Array{Int64,1} = find_index(infonames_all, variables);
	stateindex_lagged::Array{Int64,1} = find_index(statenames_lagged, variables);
	stateindex_current::Array{Int64,1} = find_index(statenames_current, variables);
	stateindex_all::Array{Int64,1} = find_index(statenames_all, variables);
	outputindex_current::Array{Int64,1} = find_index(outputnames_current, variables);
	outputindex_lead::Array{Int64,1} = find_index(outputnames_lead, variables);
	outputindex_all::Array{Int64,1} = find_index(outputnames_all, variables);
	expectindex_current::Array{Int64,1} = find_index(expectnames_current, variables);
	expectindex_lead::Array{Int64,1} = find_index(expectnames_lead, variables);
	expectindex_all::Array{Int64,1} = find_index(expectnames_all, variables);

	# Add the "lag" to the statenames_all again, to be used in the step! function
	statenames_lagged = Symbol.(string.(statenames_lagged) .* "_lag")
	statenames_all = cat(statenames_lagged, statenames_current, dims = 1)

	# indices for endog and aux needn't be split as they are all current
	endogindex::Array{Int64,1} = find_index(endognames, variables);
	auxindex::Array{Int64,1} = find_index(auxnames, variables);

	@assert (expectindex_all == cat(expectindex_current, expectindex_lead, dims = 1))

	extracted_indices = EcoNNetIndices(constant = constant,
		infoindex_lagged = infoindex_lagged, infonames_lagged = infonames_lagged,
		infoindex_current = infoindex_current, infonames_current = infonames_current,
		outputindex_current = outputindex_current, outputnames_current = outputnames_current,
		outputindex_lead = outputindex_lead, outputnames_lead = outputnames_lead,
		expectindex_current = expectindex_current, expectnames_current = expectnames_current,
		expectindex_lead = expectindex_lead, expectnames_lead = expectnames_lead,
		expectindex_all = expectindex_all, expectnames_all = expectnames_all,
		stateindex_lagged = stateindex_lagged, statenames_lagged = statenames_lagged,
		stateindex_current = stateindex_current, statenames_current = statenames_current,
		stateindex_all = stateindex_all, statenames_all = statenames_all,
		endogindex = endogindex, endognames = endognames,
		auxindex = auxindex, auxnames = auxnames);

    return extracted_indices

end

"""
A function to extract inputs for prediction
"""
function extract_inputs(s::DataFrame, tt_range::Int64, indices::EcoNNetIndices, options::EcoNNetOptions)
	# Inputs are lagged or current values of info set
	inputs::Array{Float64,1} = Array{Float64,1}(undef,0)
	if len(indices.infoindex_lagged) > 0 && len(indices.infoindex_current) > 0
		inputs = cat(Vector(s[(tt_range-1), indices.infoindex_lagged]),
			Vector(s[(tt_range), indices.infoindex_current]), dims = 1)
	elseif len(indices.infoindex_lagged) > 0
		inputs = Vector(s[(tt_range-1), indices.stateindex_lagged]);
	elseif len(indices.infoindex_current) > 0
		inputs = Vector(s[(tt_range), indices.stateindex_current]);
	else
		display(println("Need to specify some states"))
	end

	if indices.constant == 1
		inputs = cat(1.0, inputs, dims = 1);
	end

	return inputs
end
function extract_inputs(s::DataFrame, tt_range::Array{Int64,1}, indices::EcoNNetIndices)
	# Inputs are lagged or current values from information set
	inputs::Array{Float64,2} = Array{Float64,2}(undef,0,0)
	if len(indices.infoindex_lagged) > 0 && len(indices.infoindex_current) > 0
		inputs = hcat(Matrix(s[(tt_range.-1), indices.infoindex_lagged]),
			Matrix(s[(tt_range), indices.infoindex_current]));
	elseif len(indices.infoindex_lagged) > 0
		inputs = Matrix(s[(tt_range.-1), indices.infoindex_lagged]);
	elseif len(indices.infoindex_current) > 0
		inputs = Matrix(s[(tt_range), indices.infoindex_current])
	else
		display(println("Need to specify some states"))
	end
	if indices.constant == 1
		inputs = hcat(ones(size(inputs)[1]), inputs)
	end

	return inputs
end

"""
A function to extract states for step!
"""
function extract_states(s::DataFrame, tt::Int64, indices::EcoNNetIndices)
	# Inputs are lagged or current values of info set
	states::Array{Float64,1} = Array{Float64,1}(undef,0)
	if len(indices.stateindex_lagged) > 0 && len(indices.stateindex_current) > 0
		states = cat(Vector(s[(tt-1), indices.stateindex_lagged]),
			Vector(s[(tt), indices.stateindex_current]), dims = 1);
	elseif len(indices.stateindex_lagged) > 0
		states = Vector(s[(tt-1), indices.stateindex_lagged]);
	elseif len(indices.stateindex_current) > 0
		states = Vector(s[(tt), indices.stateindex_current]);
	else
		display(println("Need to specify some states"))
	end
	return states
end
function extract_states(s::DataFrame, tt::Array{Int64,1}, indices::EcoNNetIndices)
	# Inputs are lagged or current values of info set
	if len(indices.stateindex_lagged) > 0 && len(indices.stateindex_current) > 0
		states = cat(Vector(s[(tt-1), indices.stateindex_lagged]),
			Vector(s[(tt), indices.stateindex_current]), dims = 1);
	elseif len(indices.stateindex_lagged) > 0
		states = Vector(s[(tt-1), indices.stateindex_lagged]);
	elseif len(indices.stateindex_current) > 0
		states = Vector(s[(tt), indices.stateindex_current]);
	else
		display(println("Need to specify some states"))
	end
	return states
end


"""
A function to populate s with a vector (as there seems to be a bug in DataFrames)
"""
function populate(srow::DataFrameRow, index::Array{Int64,1}, newvals::Array{Float64,1})
	for cc in 1:len(index)
		srow[index[cc]] = (newvals[cc]);
	end
	return srow
end

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

#msupport, mtransition = Rouwenhorst(p,q,m,12)


"""
A function that solves a quadratic equation
"""
function solvequadratic(a, b, c)
    d = sqrt(Complex(b^2 - 4a*c))
    solns = (-b - d) / 2a, (-b + d) / 2a
    if d != real.(d)
        println("No real solutions")
    else
        solns = real.(solns)
    end
	if a == 0.0
		solns = -c/b
	elseif c == 0.0
		solns = -b/a
	end
    return solns

end




"""
Function that finds the nearest value to x in an array
	A is an array of options from which the nearest will be returned
	x is the array which needs to be matched to a row of A
"""
function findnearest(A::Array{Float64,1},x::Float64)

	indx = findmin(abs.(A.-x))[2]
	soln = A[indx]

	return soln
end



"""
Define the beliefs, predict, loss and optim function.
The two options supported here are:
(1) A feedforward neural network, provided by the Flux.jl packages
(2) Recursive Least Squares, as in Evans and Honkapohja (2001)
"""
function initialise_beliefs(options::EcoNNetOptions)
	if options.NNet
	    """
	    beliefs defined as a neural network which generates expectations.
	    """
		if options.hidden_layers == 1
	    	global beliefs = Chain(
	        	Dense(len(options.infoset), options.num_nodes, options.activation, initW = options.init_weights),
	        	Dense(options.num_nodes, len(options.expectations), initW = options.init_weights))
		elseif options.hidden_layers == 2
			global beliefs = Chain(
	        	Dense(len(options.infoset), options.num_nodes, options.activation, initW = options.init_weights),
				Dense(options.num_nodes, options.num_nodes, options.activation, initW = options.init_weights),
	        	Dense(options.num_nodes, len(options.expectations), initW = options.init_weights))
		end
	    # Here we can specify some options for the neural network
	    #loss(x, y) = mse(beliefs(x), y);

	else
	    """
	    beliefs defined as a Dict containing the parameters of an OLS regression(s) which generate expectations.
	    Each row of beta corresponds to one coefficients for one expectation.
	    The learning gain will be 1/options.window
	    """
	    beliefs = Dict([("R", Float64.(I(len(options.infoset)))),
	        ("beta" , zeros(len(options.expectations), len(options.infoset)))]);

	end
	return beliefs
end

"""
Loss function for training the network/linear regression
"""
function loss(x, y)
	mse(beliefs(x), y)
end

function predict!(inputs::Array{Float64,1}, beliefs::Chain)
	outputs::Array{Float64,1} = beliefs(inputs)
	return outputs
end
function predict!(inputs::Array{Float64,1}, beliefs::Dict)
	outputs::Array{Float64,1} = beliefs["beta"]*inputs;
	return outputs
end


"""
A function to cycle through a dictionary containing steady state values and initialise s
	Has option to alternate between two steady states ever "gap" rows.
"""
function initialise_df(df::DataFrame, steadystate::Dict{Symbol,Float64}; gap = nothing, steadystate_alt = nothing)
	newdf = deepcopy(df)
	if steadystate_alt == nothing
    	for (key, value) in steadystate
			newdf[:,key] = value*ones(nrow(df))
    	end
	elseif typeof(gap) == Int64
		#display("GAP")
		@assert nrow(df)%gap == 0 "Gap needs to be factor of nrow(df)"
		idx = Array(1:gap:nrow(newdf))
		for ii in 1:Int(nrow(newdf)/gap)
			if ii%2 == 1
				for (key, value) in steadystate
		        	newdf[idx[ii]:(idx[ii]+gap-1),key] = value*ones(gap)
		    	end
			end
			if ii%2 == 0
				for (key, value) in steadystate_alt
					newdf[idx[ii]:(idx[ii]+gap-1),key] = value*ones(gap)
		    	end
			end
		end
	end

    return newdf
end


"""
A function that names the endogenous variables, states and predictions in the step! function
"""
function name_variables(indices,x, states, predictions)
	# Name elements of x according to the endogenous variables for clarity
	# First the starting values for the endogenous variables
	for ii in 1:len(indices.endognames)
		varname = indices.endognames[ii]; value = x[ii];
		@eval (global ($varname) = ($value))
	end
	# Second the state variables
	for ii in 1:len(indices.statenames_all)
		varname = indices.statenames_all[ii]; value = states[ii];
		@eval (global ($varname) = ($value))
	end
	# Third the predictions/expectations
	for ii in 1:len(indices.expectnames_all)
		varname = indices.expectnames_all[ii]; value = predictions[ii];
		@eval (global ($varname) = ($value))
	end
end


"""
The step function determines endogeneous variables given state variables and expectations.
This function uses the equilibrium_conditions and name_variables, so is
	more convenient to use, but considerably slower due to the global variables
"""
function step!(x)
    """
    starting_values gives initial guesses for the endogeneous state variables.
    states are the state variables, both endogenous and exogenous.
    predictions is a vector containing the expectations.

    F is the output, which should be a vector of zeros
    """
	# Extract starting_values, states and predictions from input vector
	starting_values = x[1:len(options.endogenous)]
	states = x[len(options.endogenous)+1:len(options.endogenous)+len(options.states)]
	predictions = x[len(x)+1-len(options.expectations):len(x)]

    # Solve system of non-linear equations
	F = zeros(len(options.endogenous));
    results = nlsolve((F,x) -> equilibrium_conditions(F,x,states,predictions), starting_values)

    # Reduce error tolerance if solution hasn't converged
    if !(results.f_converged)
        display(println("Not converged"))
        results = nlsolve(step!, starting_values, ftol = 1.0e-12)
        display(println("New convergence: ", results.f_converged))
    end

    @assert !(NaN in Set(results.zero))

    return results.zero
end


"""
A faster version of the step function that doesn't use the name_variables function
	and so requires equilibrium conditions that manually extracts the input vector
"""
function step_fast!(x::Array{Float64,1}, options::EcoNNetOptions)
    """
    starting_values gives initial guesses for the endogeneous state variables.
    states are the state variables, both endogenous and exogenous.
    predictions is a vector containing the expectations.

    F is the output, which should be a vector of zeros
    """
	# Extract starting_values, states and predictions from input vector
	starting_values::Array{Float64,1} = x[1:len(options.endogenous)]
	states::Array{Float64,1} = x[len(options.endogenous)+1:len(options.endogenous)+len(options.states)]
	predictions::Array{Float64,1} = x[len(x)+1-len(options.expectations):len(x)]

    # Solve system of non-linear equations
	F::Array{Float64,1} = zeros(len(options.endogenous))
    results = nlsolve((F,x) -> equilibrium_conditions_fast(F,x,states,predictions), starting_values)::SolverResults
	@assert results.f_converged "Equilibrium conditions not converged - no possible solution?"
	if abs(results.residual_norm) > 0.0000001
		display("NOT CONVERGED")
	end
    #@assert !(NaN in Set(results.zero))
	new_vals::Array{Float64,1} = results.zero

    return new_vals
end


"""
The learn! function updates the beliefs (through learning)
It supports either neural network or RLS learning
It takes previous beliefs as an input, and returns updated beliefs,
"""
function learn!(beliefs::Chain, s::DataFrame, tt::Int64, options::EcoNNetOptions,
	indices::EcoNNetIndices, loss; n_cycles::Int64 = 10 )

	#println("Training network at simulation ", tt)
	# Inputs are lagged or current values of info set
	tt_range::Array{Int64,1} = Array{Int64,1}((tt-options.window+1):(tt-1))
    inputs::Array{Float64,2} = Matrix(transpose(extract_inputs(s,tt_range,indices)))
	outputs::Array{Float64,2} = Array{Float64,2}(undef,0,0)
    # Output can be both current and future values of variables
	if len(indices.outputindex_current) > 0
		outputs = hcat(Matrix(s[tt_range, indices.outputindex_current]),
        	Matrix(s[(tt_range.+1), indices.outputindex_lead]))::Array{Float64,2}
	else
		outputs = Matrix(s[(tt_range.+1), indices.outputindex_lead])::Array{Float64,2}
	end

	outputs = Matrix(transpose(outputs))

    @assert size(inputs)[2] == size(outputs)[2] "Different number of obs in inputs and outputs"

    # Split the data into training, test and validation sets
	df::NCycle = ncycle([(inputs, outputs)], n_cycles)
	if options.train_split[1] < 1.0
		props::Array{Float64,1} = [options.train_split[1], sum(options.train_split[1:2])]::Array{Float64,1}
    	idx::Array{Int64,1} = shuffle(1:size(inputs)[2])
    	train_idx::Array{Int64,1} = Array{Int64,1}(view(idx, 1:floor(Int, props[1]*size(inputs)[2])))
    	val_idx::Array{Int64,1} = Array{Int64,1}(view(idx, (floor(Int, props[1]*size(inputs)[2])+1):(floor(Int, props[2]*size(inputs)[2]))))
		test_idx::Array{Int64,1} = Array{Int64,1}(view(idx, (floor(Int, props[2]*size(inputs)[2])+1):size(inputs)[2]))

		train_inputs::Array{Float64,2} = inputs[:,train_idx]
		train_outputs::Array{Float64,2} = outputs[:,train_idx]
		val_inputs::Array{Float64,2} = inputs[:,val_idx];
		val_outputs::Array{Float64,2} = outputs[:,val_idx]
		test_inputs::Array{Float64,2} = inputs[:,test_idx];
		test_outputs::Array{Float64,2} = outputs[:,test_idx]

    	df = ncycle([(train_inputs, train_outputs)], 10)

	end

	val_loss::Array{Float64,1} = zeros(options.max_iter)
	current_loss::Float64 = 0.0
	early_stop::Int64 = 0
	prog = Progress(options.max_iter, dt = 1, desc="Progress: ")
	for ii in 1:options.max_iter

		Flux.train!(loss, params(beliefs), df, options.optim())

		if options.train_split[1] < 1.0
			current_loss = loss(val_inputs, val_outputs)
			val_loss[ii] = loss(val_inputs, val_outputs)::Float64
			#loss(train_inputs, train_outputs)
		else
			current_loss = loss(inputs, outputs)
			val_loss[ii] = current_loss
		end
		if ii > 1 && val_loss[ii] >= val_loss[ii-1]
			early_stop += 1
		else
			early_stop = 0
		end

		# Update progress bar
		next!(prog)

		if val_loss[ii] == 0.0 || early_stop >= 1
    		break
		end
	end
	display(join(["Loss: ", current_loss]))

    return beliefs
end

"""
The learn! function updates the beliefs (through learning)
It supports either neural network or RLS learning
It takes previous beliefs as an input, and returns updated beliefs,
"""
function learn!(beliefs::Dict, s, tt, options, indices, loss; optim = nothing, ss_reminder = nothing)
	### This is the RLS version
	# Note that for RLS we use one observation for each update
    # Inputs are lagged or current values of info set
    inputs = extract_inputs(s,tt-1,indices,options);
    # Outputs are current values first and leads second
	if len(indices.outputindex_current) > 0
		outputs = hcat(Vector(s[(tt-1), indices.outputindex_current]),
        	Vector(s[tt, indices.outputindex_lead]));
	else
		outputs = Vector(s[(tt), indices.outputindex_lead]);
	end
    #outputs = cat(Vector(s[(tt-1), indices.outputindex_current]),
    #    Vector(s[(tt), indices.outputindex_lead]), dims = 1);
    beliefs["R"] = beliefs["R"] + (1/options.window)*
        (inputs*transpose(inputs) - beliefs["R"])

    # Loop over each variable that needs to be predicted
    temp = (beliefs["beta"]);
    inv_R = inv(beliefs["R"]) # invert R to save multiple inversions
    for ii in 1:len(options.expectations)
        temp[ii,:] = temp[ii,:] + (1/options.window)*inv_R *inputs*
            (outputs[ii] - transpose(inputs)*temp[ii,:]);
    end
    # beta has row for each expectation, current first and lead second
    beliefs["beta"] = temp;

    return beliefs
end


"""
Function that simulates neural network learning and returns simulated data and updated beliefs
"""
function simulate_learning(sim_range::UnitRange{Int64}, s::DataFrame, beliefs::Chain, indices::EcoNNetIndices, options::EcoNNetOptions)

    # Placeholder for true expectations
    true_preds::Array{Float64,2} = ones(nrow(s),length(options.expectations))

    # Loop over the time periods specified in sim_range
    for tt in sim_range
        # Form expectations
        predictions::Array{Float64,1} = zeros(len(indices.expectnames_all))
        if tt < options.burnin && !options.burnin_use_net
            predictions = cat(Array(s[tt-1, indices.outputindex_current]), Array(s[tt-1,indices.outputindex_lead]), dims = 1);
        else
            inputs::Array{Float64,1} = extract_inputs(s,tt,indices,options);
            predictions = predict!(inputs, beliefs);
            # For numerical stability at beginning of simulation, bound the expectations
            true_preds[tt,:] = predictions
            predictions = max.(predictions, -10)
            predictions = min.(predictions, 10)
        end

        # Extract states to feed into step!
        states::Array{Float64,1} = extract_states(s,tt,indices)

        # Solve for the new values of the endogneous variables
        starting_values::Array{Float64,1} = Vector(s[tt-1,indices.endogindex]);
        new_values::Array{Float64,1} = step_fast!(cat(starting_values, states, predictions, dims = 1), options)

        # Populate the s dataframe with the new values
        s[tt,:] = populate(s[tt,:], indices.endogindex, new_values)::DataFrameRow
        s[tt,:] = populate(s[tt,:], indices.expectindex_all, predictions)::DataFrameRow

        # Update beliefs
        if tt%options.learning_gap == 0 && tt >= options.burnin
            display(join(["Simulation ", tt]))
            beliefs::Chain = learn!(beliefs, s, tt, options, indices, loss)
            #global beliefs = learn!(beliefs, s, tt, options, indices, loss, ss_reminder = DNWR)
        end

        # Plot the evolution of variables
        if options.show_plots
            if tt%options.plotting_gap == 0 && tt > (options.plotting_window+1)
                plot(legend = :bottomleft, xguidefontsize=8)
                for vv in options.plot_vars
                    (plot!(s[(tt - options.plotting_window):tt,vv], label = string(vv), legend = :bottomleft,xguidefontsize=8))
                end
                display(plot!(title = ""))
            end
        end
    end

    return beliefs,s
end





"""
Function that creates a DataFrame from a grid of initial states and beliefs
    The state_grid needs to be an N by K array where there are K state variables
"""
function create_df_grid(state_grid::Array{Float64,2}, beliefs::Chain, indices::EcoNNetIndices;
	exp_lb::Float64 = NaN, exp_ub::Float64 = NaN)
    # Create predictions from the beliefs and state_grid
    prediction_grid::Array{Float64,2} = beliefs(Matrix(transpose(state_grid)))
    prediction_grid = Matrix((transpose(prediction_grid)))
	if !isnan(exp_lb)
		prediction_grid = max.(prediction_grid, exp_lb)
	end
	if !isnan(exp_ub)
		prediction_grid = min.(prediction_grid, exp_ub)
	end

    # Create variable names
    endog_lead::Array{Symbol,1} = Symbol.(String.(indices.endognames) .* "_lead")
    exog_lead::Array{Symbol,1} = Symbol.(String.(indices.statenames_current) .* "_lead")
    grid_variables::Array{Symbol,1} = Symbol.(cat(Vector(indices.statenames_all),
        Vector(indices.expectnames_all), Vector(indices.endognames),
        endog_lead, exog_lead, dims = 1));

    # Create df_grid as DataFrame to keep strack of solution
    df_grid::DataFrame = DataFrame(ones(size(state_grid,1), length(grid_variables)), grid_variables);

    # Popoluate the first two sections of df_grid with initial state variables and predictions
    df_grid[:,indices.statenames_all] = state_grid
    df_grid[:,indices.expectnames_all] = prediction_grid

    return df_grid

end



"""
A function that takes in a df of initial conditions and runs the model
	forward two steps, taking beliefs as given.
	Additional a shock_density variable is generated which gives the probability of
	The shock which generated the t+1 variables and can be used to weight the loss function
"""
function step_map_old(df_grid::DataFrame, beliefs::Chain, indices::EcoNNetIndices, options::EcoNNetOptions,
	shock_range::Array{Float64}, shock_density::Array{Float64};
	shock_persistence::Float64 = 0.0, use_pmap::Bool = false, exp_lb::Float64 = NaN, exp_ub::Float64 = NaN)

    """
    First step: compute the period t values of endogenous variables
    """
    # Get the predictions and some starting values for nlsolve
    state_grid::Array{Float64,2} = Matrix(df_grid[:,indices.statenames_all])
    prediction_grid::Array{Float64,2} = Matrix(df_grid[:,indices.expectnames_all])

    # Use contemporaneous predictions to initialise solution for endogenous variables
	startval_grid::Array{Float64,2} = ones(size(df_grid[:,indices.endognames]))
	if len(indices.expectnames_current) == len(indices.endognames)
        startval_grid = Matrix(df_grid[:,indices.expectnames_current])
    end

    grid1::Array{Float64,2} = hcat(startval_grid, state_grid, prediction_grid)
	grid_map::Array{Array{Float64,1},1} = [grid1[i, :] for i in 1:size(df_grid,1)]
	out_grid_map::Array{Array{Float64,1},1}  = fill(Float64[],size(df_grid,1))
	# Move one step forward
    display("First step: computing period t endogenous variables")
    if !use_pmap
        @time begin
			out_grid_map = map(x -> step_fast!(x,options), grid_map)
        end
    else
        @time begin
			out_grid_map = pmap(x -> step_fast!(x,options), grid_map)
        end
    end
    out_grid::Array{Float64,2} = permutedims(reshape(hcat(out_grid_map...), (length(out_grid_map[1]), length(out_grid_map))))

    df_grid[:,indices.endognames] = out_grid

    """
    Second step: use the new states to compute the period t + 1 values
    """
    # Identify the columns that still need to be filled
    endog_lead::Array{Symbol,1} = Symbol.(String.(indices.endognames) .* "_lead")
    exog_lead::Array{Symbol,1} = Symbol.(String.(indices.statenames_current) .* "_lead")

    # The current state variables are complicated as they follow exogneous processes
    df_grid_new::DataFrame = repeat(df_grid, inner = len(shock_range))
    new_exo_mean = shock_persistence.*Matrix(df_grid_new[:,indices.statenames_current])
    df_grid_new[:,exog_lead] = new_exo_mean + repeat(shock_range,nrow(df_grid))
	#df_grid_new[:,exog_lead] = hcat(repeat(shock_range,nrow(df_grid)))

    # Different values for the new shocks have a density
    df_grid_new.shock_density = repeat(shock_density,nrow(df_grid))::Array{Float64,1}
	#df_grid_new.shock_density = shock_density

	# Remove any observations to which we assign zero probability
	df_grid_new = df_grid_new[(df_grid_new.shock_density .> 0.0),:]

	# Store realisations of contemporaneous state variables for t+1
	current_lead_grid::Array{Float64,2} = Matrix(df_grid_new[:,exog_lead])

    """
    Edit here to get the product of densities if there are multiple shocks
    """

    # Generate a fresh input grid to compute the period t+1 variables
    new_state_names::Array{Symbol,1} = Symbol.(replace.(String.(indices.statenames_lagged), "_lag" => ""))
	if len(new_state_names) > 0
		lag_lead_grid::Array{Float64,2} = Matrix(df_grid_new[:,new_state_names])
	end


    # Create new state grid to feed in for period t+1
    n_obs::Int64 = nrow(df_grid_new)
    n_states::Int64 = len(indices.statenames_all)
    state_grid_new::Array{Float64,2} = zeros((n_obs, len(indices.statenames_all)))
	if len(new_state_names) > 0
    	state_grid_new[:,1:len(indices.statenames_lagged)] = lag_lead_grid
	end
    state_grid_new[:,(len(indices.statenames_lagged)+1):n_states] = current_lead_grid

    # Create new predictions (no need to store this in the df)
    prediction_grid_new::Array{Float64,2} = beliefs(Matrix(transpose(state_grid_new)))
    prediction_grid_new = Matrix(transpose(prediction_grid_new))
	if !isnan(exp_lb)
		prediction_grid_new = max.(prediction_grid_new, exp_lb)
	end

    # Use contemporaneous predictions to initialise solution for endogenous variables
	startval_grid_new::Array{Float64,2} = ones(size(df_grid_new[:,indices.endognames]))
	if len(indices.expectnames_current) == len(indices.endognames)
        startval_grid_new = prediction_grid_new[:,1:len(indices.expectnames_current)]
    end

    # Create the new grid to pass to the map
    grid_new::Array{Float64,2} = hcat(startval_grid_new, state_grid_new, prediction_grid_new)
    grid_new_map::Array{Array{Float64,1},1} = [grid_new[i, :] for i in 1:size(grid_new,1)]
	out_grid_new_map::Array{Array{Float64,1},1}  = fill(Float64[],size(df_grid_new,1))

    # Move one step forward
    display("Second step: computing period t + 1 endogenous variables")
    if !use_pmap
        @time begin
			out_grid_new_map = map(x -> step_fast!(x,options), grid_new_map[1:100000,:])
        end
		#@time begin
		#	for ii in 1:100000
		#		out_grid_new_map[ii] = step_fast!(grid_new_map[ii],options)
		#	end
        #end
		#@time begin
		#	@sync @distributed for ii in 1:100000
		#		out_grid_new_map[ii] = step_fast!(grid_new_map[ii],options)
		#	end
        #end
    else
        @time begin
			out_grid_new_map = pmap(x -> step_fast!(x,options), grid_new)
        end
    end
    out_grid_new::Array{Float64,2} = permutedims(reshape(hcat(out_grid_new_map...), (length(out_grid_new_map[1]), length(out_grid_new_map))))

    # Insert the t+1 values into the new df
    df_grid_new[:,endog_lead] = out_grid_new



    return df_grid_new

end


"""
Map forward one step
	sval is either "E" or a Float and gives the starting value for the non-linear solver
"""
function step_1_map(df_grid::DataFrame, beliefs::Chain, indices::EcoNNetIndices,
	options::EcoNNetOptions;
	use_pmap::Bool = false, sval = "E")

    """
    First step: compute the period t values of endogenous variables
    """
    # Get the predictions and some starting values for nlsolve
    state_grid = Matrix(df_grid[:,indices.statenames_all])
    prediction_grid = Matrix(df_grid[:,indices.expectnames_all])

    # Use contemporaneous predictions to initialise solution for endogenous variables
	# Use contemporaneous predictions to initialise solution for endogenous variables
	startval_grid = ones(size(df_grid[:,indices.endognames]))
	if (len(indices.expectnames_current) == len(indices.endognames)) & (sval == "E")
        sstartval_grid = Matrix(df_grid[:,indices.expectnames_current])
    else typeof(sval) == Float64
		startval_grid = sval.*startval_grid
	end

    grid1 = hcat(startval_grid, state_grid, prediction_grid)
	grid_map = [grid1[i, :] for i in 1:size(df_grid,1)]
	out_grid_map  = fill(Float64[],size(df_grid,1))
	# Move one step forward
    display("First step: computing period t endogenous variables")
    if !use_pmap
        @time begin
			out_grid_map = map(x -> step_fast!(x,options), grid_map)
        end
		#for ii in 1:len(grid_map)
		#	display(ii)
		#	out_grid_map[ii] = step_fast!(grid_map[ii],options)
		#end
    else
        @time begin
			out_grid_map = pmap(x -> step_fast!(x,options), grid_map)
        end
    end
    out_grid = permutedims(reshape(hcat(out_grid_map...), (length(out_grid_map[1]), length(out_grid_map))))

    df_grid[:,indices.endognames] = out_grid

	return df_grid
end

"""
Function to create a transition probability matrix for a dataframe give transition ptobabilities
	transition probabilities should be a names tuple containing a DataFrame for each shock
	The first column is period t shock, the second the t+1 shock and
	the third is the corresponding transition probability
"""
function create_trans_probs(df_grid::DataFrame, transition_probabilities::NamedTuple,
	indices::EcoNNetIndices, shock_range::Array{Float64,2};parallel::Bool = false)

	# Names of the current and future shocks in the transition DataFrames
	endog_lead::Array{Symbol,1} = Symbol.(String.(indices.endognames) .* "_lead")
	exog_lead::Array{Symbol,1} = Symbol.(String.(indices.statenames_current) .* "_lead")
	exog_current::Array{Symbol,1} = indices.statenames_current

	# Create new df_grid for each possible t+1 realisation
	df_grid_new::DataFrame = repeat(df_grid, inner = size(shock_range,1))
	df_grid_new[:,exog_lead] = hcat(repeat(shock_range,nrow(df_grid)))

	# Isolate the key ingredients to calculate the shock density
	current_shocks::DataFrame = df_grid_new[:,exog_current]
	new_shocks::DataFrame = df_grid_new[:,exog_lead]
	prob_densities::DataFrame = DataFrame(zeros(size(new_shocks)))
	rename!(prob_densities, [Symbol("prob_$pp") for pp in exog_lead])
	exog_lead_probs::Array{Symbol,1} = names(prob_densities)

	# Create DataFrame for all the probabilities
	probs_df::DataFrame = hcat(current_shocks, new_shocks, prob_densities)
	probs_df.trans_prob = zeros(nrow(probs_df))

	# Cycle through each shock and fill probs_df
	for vv in 1:len(exog_current)

		current_name::Symbol = exog_current[vv]
		lead_name::Symbol = exog_lead[vv]
		prob_name::Symbol = exog_lead_probs[vv]


		if parallel
			probs_df_temp = SharedArray{Float64}(
				Array(probs_df[:,[current_name,lead_name,prob_name]]))
			transitions = SharedArray{Float64}(
				Array(transition_probabilities[current_name]))
			p_bool = SharedArray{Int64}([1])
			interval=SharedArray{Int64}([round(size(transitions,1)/5)])
		else
			probs_df_temp =
				Array(probs_df[:,[current_name,lead_name,prob_name]])
			transitions =
				Array(transition_probabilities[current_name])
			p_bool = [0]
		end


		#display(transitions)
		# M_lead
		if p_bool[1] == 0
			prog = Progress(size(transitions,1), dt = 1,
				desc=join(["Creating transition probabilities for ",exog_current[vv], ": "]))

			for ii in 1:size(transitions,1)
				# Extract the transition probability for (t,t+1) combination
				transition_prob = transitions[ii,3]
				shock_t = transitions[ii,1]
				shock_tp1 = transitions[ii,2]

				rows = (probs_df_temp[:,1] .== shock_t).*
					(probs_df_temp[:,2] .== shock_tp1)
					probs_df_temp[rows,3] .= transition_prob
				next!(prog)
			end
		elseif p_bool[1] == 1

			display(join(["Creating transition probabilities for ",exog_current[vv], " (parallel): "]))
			@sync @distributed for ii in 1:size(transitions,1)
				if ii%interval[1] == 0
					display([join([ii," out of ", size(transitions,1)])])
				end
				# Extract the transition probability for (t,t+1) combination
				transition_prob = transitions[ii,3]
				shock_t = transitions[ii,1]
				shock_tp1 = transitions[ii,2]

				rows = (probs_df_temp[:,1] .== shock_t).*
					(probs_df_temp[:,2] .== shock_tp1)

				probs_df_temp[rows,3] .= transition_prob
			end
		end

		# Insert these values into the DataFrame
		probs_df[:,prob_name] =probs_df_temp[:,3]
	end

	probs_df.trans_prob = vec(prod(Matrix(probs_df[:,exog_lead_probs]), dims = 2))::Array{Float64,1}

	@assert (df_grid_new[:,exog_current] == probs_df[:,exog_current]) "Some problem with merging period t shock"
	@assert (df_grid_new[:,exog_lead] == probs_df[:,exog_lead]) "Some problem with merging period t+1 shock"
	df_grid_new.trans_prob = probs_df.trans_prob

	# Remove any observations to which we assign zero probability
	df_grid_new = df_grid_new[(df_grid_new.trans_prob .> 0.0),:]

	# Check that the transition probs sum to the total number of initial points
	@assert round(sum(df_grid_new.trans_prob)) == round(nrow(df_grid)) "Transition probabilities don't sum to number of initial points"

	# The WSE loss function needs a (n_output, n_obs) array of weights
	loss_weights::Array{Float64,2} = Matrix(transpose(df_grid_new.trans_prob))
	loss_weights ./=(sum(loss_weights)/size(loss_weights,2))
	loss_weights = repeat(loss_weights, len(indices.expectnames_all))


	# Return the new df for grid points and also the probs_df as a sense-check
	return df_grid_new, loss_weights, probs_df
end


function step_2_prep(df_grid::DataFrame, beliefs::Chain, indices::EcoNNetIndices, options::EcoNNetOptions,
	exp_lb::Float64 = NaN, exp_ub::Float64 = NaN; sval = "E")

    """
    Second step: use the new states to compute the period t + 1 values
    """
    # Identify the columns that still need to be filled
    endog_lead::Array{Symbol,1} = Symbol.(String.(indices.endognames) .* "_lead")
    exog_lead::Array{Symbol,1} = Symbol.(String.(indices.statenames_current) .* "_lead")

	# Store realisations of contemporaneous state variables for t+1
	current_lead_grid::Array{Float64,2} = Matrix(df_grid_new[:,exog_lead])

    # Generate a fresh input grid to compute the period t+1 variables
    new_state_names::Array{Symbol,1} = Symbol.(replace.(String.(indices.statenames_lagged), "_lag" => ""))
	if len(new_state_names) > 0
		lag_lead_grid::Array{Float64,2} = Matrix(df_grid_new[:,new_state_names])
	end

	# Create new state grid to feed in for period t+1
    n_obs::Int64 = nrow(df_grid_new)
    n_states::Int64 = len(indices.statenames_all)
    state_grid_new::Array{Float64,2} = zeros((n_obs, len(indices.statenames_all)))
	if len(new_state_names) > 0
    	state_grid_new[:,1:len(indices.statenames_lagged)] = lag_lead_grid
	end
    state_grid_new[:,(len(indices.statenames_lagged)+1):n_states] = current_lead_grid

    # Create new predictions (no need to store this in the df)
    prediction_grid_new::Array{Float64,2} = beliefs(Matrix(transpose(state_grid_new)))
    prediction_grid_new = Matrix(transpose(prediction_grid_new))
	if !isnan(exp_lb)
		prediction_grid_new = max.(prediction_grid_new, exp_lb)
	end
	if !isnan(exp_ub)
		prediction_grid_new = min.(prediction_grid_new, exp_ub)
	end

    # Use contemporaneous predictions to initialise solution for endogenous variables
	startval_grid_new = ones(size(df_grid_new[:,indices.endognames]))
	if (len(indices.expectnames_current) == len(indices.endognames)) & (sval == "E")
        startval_grid_new = prediction_grid_new[:,1:len(indices.expectnames_current)]
    else typeof(sval) == Float64
		startval_grid_new = sval.*startval_grid_new
	end
    # Create the new grid to pass to the map
    grid_new::Array{Float64,2} = hcat(startval_grid_new, state_grid_new, prediction_grid_new)
    #grid_new_map::Array{Array{Float64,1},1} = [grid_new[i, :] for i in 1:size(grid_new,1)]

	return grid_new, endog_lead

end



"""
Function that defines weighted mean square error loss for beliefs
	Need to specify weights::Array{Float64,2}(size(outputs))
"""
function WMSE(x::Array{Float64,2}, y::Array{Float64,2}, weights::Array{Float64,2})
	se = (beliefs(x) .- y).^2
	wse = weights.*se
	wmse_loss = sum(wse) * 1 // length(y)
	return wmse_loss
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
		CP[trans_probs[st,1],trans_probs[st,2]] = trans_probs[st,3]
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



"""
Function to classify state at t+1 in df_grid_new
"""
function classify_states(df_grid, df_grid_new,range_dicts,state_tp1_vars)
	range_lengths = []
	for vv in 1:len(range_dicts)
		push!(range_lengths,len(range_dicts[vv]))
	end
	prog = Progress(nrow(df_grid_new), dt = 1, desc="Classifying t+1 state: ")
	for df_row in 1:nrow(df_grid_new)
		#state_defs = unique(df_grid_new[:,state_t_vars])
		state_def = df_grid_new[df_row,state_tp1_vars]

		# Find the index of each variable
		π_def = π_dict[state_def.π]
		y_def = y_dict[state_def.y]
		ϵ_y_def = ϵ_y_dict[state_def.ϵ_y_lead]

		st = 0

		st += (ϵ_y_def-1)*len(π_range)*len(y_range)
		st += (y_def-1)*len(π_range)
		st += π_def

		@assert state_labels[st] == Array(state_def)

		# Label state at t+1
		df_grid_new[df_row,:state_tp1]=st
		# Use this to evaluate endogenous variables
		endog_new = df_grid[st,endog]
		for vv in 1:len(endog_new)
			df_grid_new[df_row,endog_lead[vv]]= endog_new[vv]
		end
		next!(prog)
	end

	return df_grid_new

end







"""
Function that updates the beliefs based on df_grid_new and weighted loss function
"""
function update_beliefs(df::DataFrame, beliefs::Chain, indices::EcoNNetIndices, options::EcoNNetOptions;
	epochs::Int64 = 1, cycles::Int64 = 10, verbose::Bool = true, weights::Array{Float64,2} = ones(2,2))

	# Extract inputs and outputs
    outputnames_new = Symbol.(replace.(String.(indices.expectnames_all), "E"=> ""))
    outputs = Matrix(transpose(Matrix(df[:,outputnames_new])))
    inputs::Array{Float64,2} = Matrix(transpose(Matrix(df[:,indices.statenames_all])))

	idx::Array{Int64,1} = shuffle(1:size(inputs)[2])

	inputs = inputs[:,idx]
	outputs = outputs[:,idx]

	train_df::NCycle = ncycle([(inputs, outputs)], cycles)

	# Update weights if none are specified
	if weights == ones(2,2)
		weights = ones(size(outputs))
	else
		weights =weights[:,idx]
	end

	initial_loss = WMSE(inputs,outputs,weights)::Float64
	display(join(["Initial loss: ", initial_loss]))

	loss_iters = zeros(epochs)
	for ii in 1:epochs
    	Flux.train!((x,y) -> WMSE(x,y,weights), params(beliefs), train_df, options.optim());
		loss_iters[ii] = WMSE(inputs,outputs,weights)
		if verbose
			display(join(["Iter: ", ii, ", Loss: ", loss_iters[ii]]))
		end
		if ii > 1
			if (loss_iters[ii] > loss_iters[ii-1])
				break
			end
		end
	end
    display(join(["Weighted forecast error is ", WMSE(inputs, outputs, weights)]))

    return beliefs
end

















"""
Function that generates a perfect foresight path from given starting point
"""
function pf_path(initial_ss; periods = 100)

	global paths = DataFrame(zeros(periods, len(variables)), variables);
    paths = initialise_df(paths, initial_ss);

    for tt in 3:periods
        #println(tt)
        inputs = vcat(Array(paths[tt-2,indices.endognames]),Array(paths[tt-1,indices.endognames]));
        new_values = perfect_foresight(inputs)
        # Populate the paths DataFrame
		paths[tt,:] = populate(paths[tt,:], indices.endogindex, new_values)::DataFrameRow
    end

    return paths

end


"""
A function that uses the beliefs and equilibrium conditions to generate IRFs
"""
function irf(type, initial_ss, beliefs::Chain; periods = 100, magnitude = 1.0, persistence = 0.9,
		show_plot = true, plot_vars = nothing, shock_period = 3, y_lim::Array{Float64,1} = zeros(2))
    global paths = DataFrame(zeros(periods, len(variables)), variables);
    paths = initialise_df(paths, initial_ss);
    global shock = zeros(periods)
    shock[shock_period] = magnitude
    shock = simulate_ar(persistence, 0.00, periods, shock)
    @eval (global (paths.$type = shock))
    #plot(paths.ϵ_y)

    for tt in shock_period:periods
        #println(tt)
        # inputs, prediction and state
        inputs = extract_inputs(paths,tt,indices,options);
        predictions = predict!(inputs, beliefs);
        states = extract_states(paths,tt,indices)

        # Solve system of non-linear equations
        starting_values = Vector(paths[tt-1,indices.endogindex]);
        new_values = step_fast!(cat(starting_values, states, predictions, dims = 1), options)

        # Populate the paths DataFrame
		paths[tt,:] = populate(paths[tt,:], indices.endogindex, new_values)::DataFrameRow
        paths[tt,:] = populate(paths[tt,:], indices.expectindex_all, predictions)::DataFrameRow

    end
	if show_plot
		#plotly()
		if y_lim == zeros(2)
    		plt = plot(paths[:,type], label = String(type))
		else
			plt = plot(paths[:,type], label = String(type),  ylims =y_lim)
		end
		for vv in 1:len(plot_vars)
			if vv < len(plot_vars)
				plot!(paths[:,plot_vars[vv]], label = plot_vars[vv])
			else
				display(plot!(paths[:,plot_vars[vv]], label = plot_vars[vv]))
			end
		end
	end

    return paths

end

"""
Plot arrows on phase diagram (at arrow_point to arrow_point+1)
"""
function phase_arrow_plot(paths, vars; arrow_points=[0],
	h_points = 11:99, v_points = 12:100)

	h_path = paths[h_points,vars[1]]
	v_path =paths[v_points,vars[2]]

	plt = plot!(h_path,v_path,color = :blue, label = false)
	for point in arrow_points
		plot!([h_path[point],h_path[point+1]], [v_path[point],v_path[point+1]],
			arrow = :closed, color = :blue, label = false)
	end
	display(plot!())
	return plt
end
