# Import dependencies
using Distributions, TableView, GLM
using Flux, Distributed, JLD
using Flux: throttle, params, mse, glorot_uniform, @epochs
using IterTools: ncycle
using DataFrames, Random, Atom
using Plots, NLsolve

"""
General convention on ordering variables:
	1) Order in chronological order (i.e. lagged, current, lead)
	2) Order alphabetically within each period
"""


"""
A shortened name for the length function to make code more compact
"""
function len(x)
	length(x)
end

"""
A function that identifies the position of elements in an array while preserving their order
"""
function find_index(element_names, variables)
	element_index = fill(0, len(element_names))
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
function infoindices(options, variables)
	infonames_lagged = Any[];
	infonames_current = Any[];
	infonames_all = Any[];
	outputnames_current = Any[];
	outputnames_lead = Any[];
	outputnames_all = Any[]
	statenames_lagged = Any[];
	statenames_current = Any[];
	statenames_all = Any[];
	expectnames_current = Any[];
	expectnames_lead = Any[];
	expectnames_all = Any[];

	constant = 0;

	auxnames = options.auxiliary;
	endognames = options.endogenous;

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

	infoindex_lagged = find_index(infonames_lagged, variables);
	infoindex_current = find_index(infonames_current, variables);
	infoindex_all = find_index(infonames_all, variables);
	stateindex_lagged = find_index(statenames_lagged, variables);
	stateindex_current = find_index(statenames_current, variables);
	stateindex_all = find_index(statenames_all, variables);
	outputindex_current = find_index(outputnames_current, variables);
	outputindex_lead = find_index(outputnames_lead, variables);
	outputindex_all = find_index(outputnames_all, variables);
	expectindex_current = find_index(expectnames_current, variables);
	expectindex_lead = find_index(expectnames_lead, variables);
	expectindex_all = find_index(expectnames_all, variables);

	# Add the "lag" to the statenames_all again, to be used in the step! function
	statenames_lagged = Symbol.(string.(statenames_lagged) .* "_lag")
	statenames_all = cat(statenames_lagged, statenames_current, dims = 1)

	# indices for endog and aux needn't be split as they are all current
	endogindex = find_index(endognames, variables);
	auxindex = find_index(auxnames, variables);

	@assert (expectindex_all == cat(expectindex_current, expectindex_lead, dims = 1))

    return (constant = constant,
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

end

"""
A function to extract inputs for prediction
"""
function extract_inputs(s, tt_range, indices, options)
	# Inputs are lagged or current values of info set
	if len(tt_range) == 1
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
	elseif len(tt_range) > 1
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
	else
		display("Need to specify valid time indices")
	end

	return inputs
end

"""
A function to extract states for step!
"""
function extract_states(s, tt, indices)
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
function populate(s, tt, index, newvals)
	if len(tt) == 1
		for cc in 1:len(index)
			global s[tt,index[cc]] = (newvals[cc]);
		end
	end
end

"""
A function that generates n draws from an AR process with persistence
	ρ and noise parameterised by σ drawn from dist distribution.
"""
function simulate_ar(ρ, σ, n, dist)
	if len(dist) == n
		noise = dist
	else
		noise = σ*rand(dist, n)
	end
	y = [0.0 for ii = 1:n]
	y[1] = noise[1]
	for ii in 2:(n)
		y[ii] = ρ*y[ii-1] + noise[ii]
	end
	return y
end;

"""
A function that creates rare large shocks, of alternating sign.
    gap is the gap between shocks
    mag is the magnitude of the shock (in absolute terms)
    start is the first period in which the shock hits
    N is the total length of the shock series
"""
function alternating_shocks(gap, mag, start, N)
    noise = zeros(N)
    noise[start:(2*gap):N] = mag*ones(len(start:(2*gap):N))
    noise[(start+gap):(2*gap):N] = -mag*ones(len((start+gap):(2*gap):N))
    return noise
end
#dist = Normal()
#dist = Multinomial(1, [0.001, 0.998, 0.001] )
#noise_d = Matrix(Float64.(rand(dist, nrow(s))))
#noise_d[1,:] = -shock_size*noise_d[1,:]
#noise_d[2,:] = 0.0*noise_d[2,:]
#noise_d[3,:] = shock_size*noise_d[3,:]
#noise_d = vec(sum(noise_d, dims = 1))



"""
A function that returns support and transition kernel for an N-state
	Markov chain Rouwenhorst approximation of an AR(1) process
	p:
	q:
	m: must be a positive scalar
	N: number of states, must be a positive integer
"""
function Rouwenhorst(p::Float64,q::Float64,m::Float64,N::Int)

	vtest = [p, q, m, N]
	@assert (size(vtest,1) == 4) "Rouwenhorst(): argument dimensions incorrect -- must all be scalars."
	@assert (p >= 0 && p < 1) "Rouwenhorst(): p must be between zero and 1."
	@assert (q >= 0 && q < 1) "Rouwenhorst(): q must be between zero and 1."
	@assert (m > 0) "Rouwenhorst(): m must be positive."
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
	vm = -m : 2*m/(N-1) : m
	vm = Array(vm)

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
A function that finds the nearest value to x in an array
"""
function findnearest(A,x::Float64)
	if len(A) > 1
		A = collect(A)
		indx = findmin(abs.(A.-x))[2]
		soln = Float64(A[indx])
	else
		soln = A
	end
	return soln
end



"""
Define the beliefs, predict, loss and optim function.
The two options supported here are:
(1) A feedforward neural network, provided by the Flux.jl packages
(2) Recursive Least Squares, as in Evans and Honkapohja (2001)
"""
function initialise_beliefs(options)
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
	    global optim = options.optim();

	else
	    """
	    beliefs defined as a Dict containing the parameters of an OLS regression(s) which generate expectations.
	    Each row of beta corresponds to one coefficients for one expectation.
	    The learning gain will be 1/options.window
	    """
	    global beliefs = Dict([("R", Float64.(I(len(options.infoset)))),
	        ("beta" , zeros(len(options.expectations), len(options.infoset)))]);

	end
end

"""
Loss function for training the network/linear regression
"""
function loss(x, y)
	if options.NNet
		mse(beliefs(x), y)
	else
		loss(x, y) = mse(predict!(x), y);
	end
end

function predict!(inputs, beliefs)
	if options.NNet
		outputs = beliefs(inputs)
	else
		outputs = beliefs["beta"]*inputs;
	end
	return outputs
end


"""
A function to cycle through a dictionary containing steady state values and initialise s
	Has option to alternate between two steady states ever "gap" rows.
"""
function initialise_s(df, steadystate; gap = nothing, steadystate_alt = nothing)
	global tempdf = df
	global gap1 = gap
	#display(gap1)
	if steadystate_alt == nothing
    	for (key, value) in steadystate
        	colname = Symbol(key)
        	@eval ((tempdf.$colname) = ($value)*ones(nrow(tempdf)))
    	end
	elseif typeof(gap1) == Int64
		@assert nrow(tempdf)%gap1 == 0
		global jj = 1
		for ii in 1:Int(nrow(tempdf)/gap1)
			if ii%2 == 1
				for (key, value) in steadystate
		        	colname = Symbol(key)
		        	@eval ((tempdf.$colname[jj:(jj+gap1-1)])) = ($value)*ones(gap1)
		    	end
			end
			if ii%2 == 0
				for (key, value) in steadystate_alt
		        	colname = Symbol(key)
		        	@eval ((tempdf.$colname[jj:(jj+gap1-1)])) = ($value)*ones(gap1)
		    	end
			end
			global jj = jj + gap1
		end
	end

    return tempdf
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
function step_fast!(x::Array{Float64,1})
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
	F = zeros(len(options.endogenous))::Array{Float64,1}
    results = nlsolve((F,x) -> equilibrium_conditions_fast(F,x,states,predictions), starting_values)
	@assert results.f_converged "Equilibrium conditions not converged - no possible solution?"
    @assert !(NaN in Set(results.zero))

    return results.zero
end


"""
The learn! function updates the beliefs (through learning)
It supports either neural network or RLS learning
It takes previous beliefs as an input, and returns updated beliefs,
"""
function learn!(beliefs::Chain, s::DataFrame, tt::Int64, options::NamedTuple,
	indices::NamedTuple, loss; optim = nothing, ss_reminder = nothing)

	if optim == nothing
		optim = options.optim()
	end
	#println("Training network at simulation ", tt)
	# Inputs are lagged or current values of info set
	tt_range = (tt-options.window+1):(tt-1)::Int64
	inputs = extract_inputs(s,tt_range,indices, options)::Array{Float64,2}
    if ss_reminder != nothing
		println("Steady state reminder included")
		s[1:1005,:] = initialise_s(s[1:1005,:], ss_reminder)
		s[1:1005,:] = initialise_s(s[1:1005,:], ss_reminder)
		temp = hcat(Matrix(s[(1005-1000):(1005-2), indices.infoindex_lagged]),
            Matrix(s[(1005-1000+1):(1005-1), indices.infoindex_current]))
		if indices.constant == 1
			temp = hcat(ones(size(temp)[1]), temp)
	    end
		inputs = vcat(temp, inputs)
	end

    inputs = Matrix(transpose(inputs))
    # Output can be both current and future values of variables
	if len(indices.outputindex_current) > 0
		outputs = hcat(Matrix(s[(tt-options.window+1):(tt-1), indices.outputindex_current]),
        	Matrix(s[(tt-options.window+2):(tt), indices.outputindex_lead]))
	else
		outputs = Matrix(s[(tt-options.window+2):(tt), indices.outputindex_lead])
	end
	if ss_reminder != nothing
		if len(indices.outputindex_current) > 0
			temp = hcat(Matrix(s[(1005-1000+1):(1005-1), indices.outputindex_current]),
            	Matrix(s[(1005-1000+2):(1005), indices.outputindex_lead]))
		else
			temp = Matrix(s[(1005-1000+2):(1005), indices.outputindex_lead])
		end
		outputs = vcat(temp, outputs)
	end
	outputs = Matrix(transpose(outputs))

    @assert size(inputs)[2] == size(outputs)[2]

    # Split the data into training, test and validation sets
	if options.train_split[1] < 1.0
		props = [options.train_split[1], sum(options.train_split[1:2])]
    	idx = shuffle(1:size(inputs)[2])
    	train_idx = view(idx, 1:floor(Int, props[1]*size(inputs)[2]))
    	val_idx = view(idx, (floor(Int, props[1]*size(inputs)[2])+1):(floor(Int, props[2]*size(inputs)[2])))
		test_idx = view(idx, (floor(Int, props[2]*size(inputs)[2])+1):size(inputs)[2])

		train_inputs = inputs[:,train_idx]
		train_outputs = outputs[:,train_idx]
		val_inputs = inputs[:,val_idx]; val_outputs = outputs[:,val_idx]
		test_inputs = inputs[:,test_idx]; test_outputs = outputs[:,test_idx]

    	df = ncycle([(train_inputs, train_outputs)], 10)
	else
		df = ncycle([(inputs, outputs)], 10)
	end

	val_loss = zeros(options.max_iter);
	global early_stop = 0;
	for ii in 1:options.max_iter

		Flux.train!(loss, params(beliefs), df, optim);

		if options.train_split[1] < 1.0
			val_loss[ii] = loss(val_inputs, val_outputs);
			#loss(train_inputs, train_outputs)
		else
			val_loss[ii] = loss(inputs, outputs)
		end
		if ii > 1 && val_loss[ii] >= val_loss[ii-1]
			global early_stop += 1
		else
			global early_stop = 0
		end

		println(ii, "Loss: ", val_loss[ii], ", Early stop: ", early_stop)

		if val_loss[ii] == 0.0 || early_stop >= 1
    		break
		end
	end

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
Function that creates a DataFrame from a grid of initial states and beliefs
    The state_grid needs to be an N by K array where there are K state variables
"""
function create_df_grid(state_grid, beliefs, indices; exp_lb = -100.0)
    # Create predictions from the beliefs and state_grid
    prediction_grid = beliefs(Matrix(transpose(state_grid)))
    prediction_grid = Matrix((transpose(prediction_grid)))
	prediction_grid = max.(prediction_grid, exp_lb)

    # Create variable names
    endog_lead = Symbol.(String.(indices.endognames) .* "_lead")
    exog_lead = Symbol.(String.(indices.statenames_current) .* "_lead")
    grid_variables = Symbol.(cat(Vector(indices.statenames_all),
        Vector(indices.expectnames_all), Vector(indices.endognames),
        endog_lead, exog_lead, dims = 1));

    # Create df_grid as DataFrame to keep strack of solution
    df_grid = DataFrame(ones(size(state_grid,1), length(grid_variables)), grid_variables);

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
function step_map(df_grid, beliefs, indices, shock_range, shock_density;
	shock_persistence = 0.0, use_pmap = false, exp_lb = -100.0)

    """
    First step: compute the period t values of endogenous variables
    """
    # Get the predictions and some starting values for nlsolve
    state_grid = Matrix(df_grid[:,indices.statenames_all])
    prediction_grid = Matrix(df_grid[:,indices.expectnames_all])

    # Use contemporaneous predictions to initialise solution for endogenous variables
    if len(indices.expectnames_current) == len(indices.endognames)
        startval_grid = Matrix(df_grid[:,indices.expectnames_current])
    else
        startval_grid = ones(size(df_grid[:,indices.endognames]))
    end

    grid = hcat(startval_grid, state_grid, prediction_grid)
    grid = [grid[i, :] for i in 1:size(grid,1)]

    # Move one step forward
    display("First step: computing period t endogenous variables")
    if !use_pmap
        @time begin
            out_grid = map(step_fast!, grid)
        end
    else
        @time begin
            out_grid = pmap(step_fast!, grid)
        end
    end
    out_grid = permutedims(reshape(hcat(out_grid...), (length(out_grid[1]), length(out_grid))))

    df_grid[:,indices.endognames] = out_grid


    """
    Second step: use the new states to compute the period t + 1 values
    """
    # Identify the columns that still need to be filled
    endog_lead = Symbol.(String.(indices.endognames) .* "_lead")
    exog_lead = Symbol.(String.(indices.statenames_current) .* "_lead")

    # The current state variables are complicated as they follow exogneous processes
    df_grid_new = repeat(df_grid, inner = len(shock_range))
    #new_exo_mean = shock_persistence.*Matrix(df_grid_new[:,indices.statenames_current])
    #df_grid_new[:,exog_lead] = new_exo_mean + repeat(shock_range,nrow(df_grid))
	df_grid_new[:,exog_lead] = hcat(repeat(shock_range,nrow(df_grid)))
    current_lead_grid = Matrix(df_grid_new[:,exog_lead])

    # Different values for the new shocks have a density
    #df_grid_new.shock_density = repeat(shock_density,nrow(df_grid))
	df_grid_new.shock_density = shock_density
    """
    Edit here to get the product of densities if there are multiple shocks
    """

    # Generate a fresh input grid to compute the period t+1 variables
    new_state_names = Symbol.(replace.(String.(indices.statenames_lagged), "_lag" => ""))
	if len(new_state_names) > 0
		lag_lead_grid = Matrix(df_grid_new[:,new_state_names])
	end


    # Create new state grid to feed in for period t+1
    n_obs = nrow(df_grid_new)
    n_states = len(indices.statenames_all)
    state_grid_new = zeros((n_obs, len(indices.statenames_all)))
	if len(new_state_names) > 0
    	state_grid_new[:,1:len(indices.statenames_lagged)] = lag_lead_grid
	end
    state_grid_new[:,(len(indices.statenames_lagged)+1):n_states] = current_lead_grid

    # Create new predictions (no need to store this in the df)
    prediction_grid_new = beliefs(Matrix(transpose(state_grid_new)))
    prediction_grid_new = Matrix(transpose(prediction_grid_new))
	prediction_grid_new = max.(prediction_grid_new, exp_lb)

    # Use contemporaneous predictions to initialise solution for endogenous variables
    if len(indices.expectnames_current) == len(indices.endognames)
        startval_grid_new = prediction_grid_new[:,1:len(indices.expectnames_current)]
    else
        startval_grid_new = ones(size(df_grid_new[:,indices.endognames]))
    end

    # Create the new grid to pass to the map
    grid_new = hcat(startval_grid_new, state_grid_new, prediction_grid_new)
    grid_new = [grid_new[i, :] for i in 1:size(grid_new,1)]

    # Move one step forward
    display("Second step: computing period t + 1 endogenous variables")
    if !use_pmap
        @time begin
            out_grid_new = map(step_fast!, grid_new)
        end
    else
        @time begin
            out_grid_new = pmap(step_fast!, grid_new)
        end
    end
    out_grid_new = permutedims(reshape(hcat(out_grid_new...), (length(out_grid_new[1]), length(out_grid_new))))

    # Insert the t+1 values into the new df
    df_grid_new[:,endog_lead] = out_grid_new

	# Remove any observations to which we assign zero probability
	df_grid_new = df_grid_new[(df_grid_new.shock_density .> 0.0),:]

    return df_grid_new

end

"""
Function that updates the beliefs based on df_grid_new and weighted loss function
"""
function update_beliefs(df, beliefs, indices, epoch; verbose = true)
	# Extract outputs
    outputnames_new = Symbol.(replace.(String.(indices.expectnames_all), "E"=> ""))
    outputs = Matrix(df[:,outputnames_new])
    outputs = Matrix(transpose(outputs))

	# Extract inputs
    inputs = Matrix(df[:,indices.statenames_all])
    inputs = Matrix(transpose(inputs))

    train_df = ncycle([(inputs, outputs)], 10)
	if verbose
		@epochs epoch Flux.train!(loss_weighted, params(beliefs), train_df, optim);
	else
		for ii in 1:epoch
	    	Flux.train!(loss_weighted, params(beliefs), train_df, optim);
		end
	end
    display(println("Weighted forecast error is ", loss_weighted(inputs, outputs)))

    return beliefs
end


















"""
A function that uses the beliefs and equilibrium conditions to generate IRFs
"""
function irf(type, initial_ss; periods = 100, magnitude = 1.0, persistence = 0.9,
		show_plot = true, plot_vars = nothing, shock_period = 3)
    global paths = DataFrame(zeros(periods, len(variables)), variables);
    paths = initialise_s(paths, initial_ss);
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
        new_values = step_fast!(cat(starting_values, states, predictions, dims = 1))

        # Populate the paths DataFrame
        populate(paths, tt, indices.endogindex, new_values)
        populate(paths, tt, indices.expectindex_all, predictions);
    end
	if show_plot
		plotly()
    	plot(paths[:,type], label = String(type))
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
