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
Function that creates a DataFrame from a grid of initial states and beliefs
    The state_grid needs to be an N by K array where there are K state variables
"""
function create_df_grid(state_grid::Array{Float64,2}, beliefs::Chain, indices::EcoNNetIndices;
	exp_lb::Float64 = NaN, exp_ub::Float64 = NaN)
    # Create predictions from the beliefs and state_grid
    prediction_grid::Array{Float64,2} = Tracker.data(beliefs(Matrix(transpose(state_grid))))
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
Function that updates the beliefs based on df_grid_new and weighted loss function
"""
function update_beliefs(df::DataFrame, beliefs::Chain, indices::EcoNNetIndices, options::EcoNNetOptions;
	epochs::Int64 = 1, cycles::Int64 = 10, verbose::Bool = true, weights::Array{Float64,2} = ones(2,2))

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

	initial_loss::Float64 = Tracker.data(WMSE(inputs,outputs,weights))
	display(join(["Initial loss: ", initial_loss]))

	loss_iters = zeros(epochs)
	for ii in 1:epochs
    	Flux.train!((x,y) -> WMSE(x,y,weights), params(beliefs), train_df, options.optim());
		loss_iters[ii] = Tracker.data(WMSE(inputs,outputs,weights))
		if verbose
			display(join(["Iter: ", ii, ", Loss: ", loss_iters[ii]]))
		end
		if ii > 1
			if (loss_iters[ii] > loss_iters[ii-1])
				break
			end
		end
	end
    display(join(["Weighted forecast error is ", Tracker.data(WMSE(inputs, outputs, weights))]))

    return beliefs
end
