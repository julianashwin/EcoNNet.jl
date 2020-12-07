"""
Function to get state t+1 variables from state t variables
"""
function state_tp1(state_t_vars)
	state_tp1_vars = String.(state_t_vars)
	for vv in 1:len(state_t_vars)
		st_str = string(state_t_vars[vv])
		if occursin("_lag", st_str)
			st_str = replace(st_str, "_lag" => "")
		else
			st_str*="_lead"
		end
		state_tp1_vars[vv] = st_str
	end
	return Symbol.(state_tp1_vars)
end


"""
Function that maps endogenous variables to their closest grid point, given by ranges
"""
function closest_gridpoint(df_grid::DataFrame, endog::Array{Symbol,1}, ranges::NamedTuple)::DataFrame
	for vv in 1:len(endog)
		df_grid[:,endog[vv]] = map(x -> findnearest(ranges[endog[vv]],x), Array(df_grid[:,endog[vv]]))
	end
	return df_grid
end


"""
Function that classifies the state at t+1 for each grid point
	The values of the endogenous variables are automatically filled, to avoid
	using the equilbrium conditions to calculate them
"""
function classify_states(df_grid_new::DataFrame, range_dicts::NamedTuple, range_lengths::NamedTuple,
		options::EcoNNetOptions)::DataFrame
	state_t_vars::Array{Symbol,1} = options.states
	state_tp1_vars::Array{Symbol,1}  = state_tp1(state_t_vars)
	class_vars::Array{Symbol,1}  = Symbol.(replace.(string.(state_tp1_vars), "_lead" => ""))
	nvars::Int64 = len(class_vars)
	defs::Array{Int64,1} = Int.(zeros(len(class_vars)))

	# Check every row of the dataframe
	prog = Progress(nrow(df_grid_new), dt = 1, desc="Classifying t+1 state: ")
	for df_row in 1:nrow(df_grid_new)

		# Find the index of each variable
		state_def = df_grid_new[df_row,state_tp1_vars]
		for vv in 1:nvars
			defs[vv] = range_dicts[class_vars[vv]][state_def[state_tp1_vars[vv]]]
		end
		# Combine these individual indices to indentify the state
		st = defs[1]
		for vv in 2:nvars
			st += (defs[vv] - 1)*prod([range_lengths[pp] for pp in class_vars[1:(vv-1)]])
		end
		#@assert state_labels[st] == Array(state_def)
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
Function that computes expectations from transition probabilities
"""
function compute_expectations(df_grid, df_grid_new, shock_range)
		start_obs = 1
		end_obs = size(shock_range,1)
		prog = Progress(nrow(df_grid), dt = 1, desc="Computing Expectations: ")
		for st in 1:nrow(df_grid)
			df_rows = df_grid_new[start_obs:end_obs,:]
			df_rows.trans_prob ./= sum(df_rows.trans_prob)
			df_grid[st,:π_lead] = sum(df_rows.trans_prob.*df_rows.π_lead)
			df_grid[st,:y_lead] = sum(df_rows.trans_prob.*df_rows.y_lead)
			start_obs += size(shock_range,1)
			end_obs += size(shock_range,1)
			next!(prog)
		end
		return df_grid
end





"""
Function that updates the beliefs based on df_grid_new and weighted loss function
"""
function update_beliefs(df::DataFrame, beliefs::Chain, indices::EcoNNetIndices, options::EcoNNetOptions;
	epochs::Int64 = 1, cycles::Int64 = 10, verbose::Bool = true, weights::Array{Float64,2} = ones(2,2),
	cutoff::Bool = true)

	# Extract inputs and outputs
    outputnames_new = Symbol.(replace.(String.(indices.expectnames_all), "E"=> ""))
    outputs::Array{Float64,2} = Matrix(transpose(Matrix(df[:,outputnames_new])))
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

	initial_loss::Float64 = WMSE(inputs,outputs,weights)
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
				if cutoff
					break
				end
			end
		end
	end
    display(join(["Weighted forecast error is ", WMSE(inputs, outputs, weights)]))

    return beliefs
end
