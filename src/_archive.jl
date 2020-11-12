
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
