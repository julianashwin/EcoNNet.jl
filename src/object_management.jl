"""
Functions to manage the various objects used throughout and extract relevant information
"""


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
