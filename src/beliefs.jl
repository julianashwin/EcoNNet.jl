"""
Functions and objects related to the beliefs neural network
"""



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
Generate prediction from beliefs and inputs
"""
function predict!(inputs::Array{Float64,1}, beliefs::Chain)
	outputs::Array{Float64,1} = beliefs(inputs)
	return outputs
end
function predict!(inputs::Array{Float64,1}, beliefs::Dict)
	outputs::Array{Float64,1} = beliefs["beta"]*inputs;
	return outputs
end





"""
The learn! function updates the beliefs (through learning)
It supports either neural network or RLS learning
It takes previous beliefs as an input, and returns updated beliefs,
"""
function learn!(beliefs::Chain, s::DataFrame, tt::Int64, options::EcoNNetOptions,
	indices::EcoNNetIndices, loss; n_cycles::Int64 = 10, cutoff::Bool = true)

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
			if cutoff
    			break
			end
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
