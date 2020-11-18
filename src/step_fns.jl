"""
Functions to solve the equilibrium conditions, given beliefs
"""


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
        begin
			out_grid_map = map(x -> step_fast!(x,options), grid_map)
        end
		#for ii in 1:len(grid_map)
		#	display(ii)
		#	out_grid_map[ii] = step_fast!(grid_map[ii],options)
		#end
    else
        begin
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
	exog_lead_probs::Array{Symbol,1} = [Symbol("prob_$pp") for pp in exog_lead]
	@assert len(exog_current) == len(transition_probabilities) "Transition probs don't match the shocks"

	# Create new df_grid for each possible t+1 realisation
	df_grid_new::DataFrame = repeat(df_grid, inner = size(shock_range,1))
	df_grid_new[:,exog_lead] = hcat(repeat(shock_range,nrow(df_grid)))
	df_grid_new.trans_prob = zeros(nrow(df_grid_new))

	# Create DataFrame for each shock's transition probability
	prob_densities::DataFrame = DataFrame(zeros(size(df_grid_new[:,exog_lead])))
	rename!(prob_densities, exog_lead_probs)

	nshocks::Int64 = len(transition_probabilities)
	range_len::Array{Int,1} = Array{Int,1}(Int.(zeros(nshocks)))
	for vv in 1:len(exog_current)
		range_len[vv] = sqrt(nrow(transition_probabilities[exog_current[vv]]))
	end

	endog_nstates::Int64 = Int(nrow(df_grid_new)/size(shock_range,1)^2)
	@assert endog_nstates*size(shock_range,1) == nrow(df_grid) "df_grid rows and nstates don't match up"

	prog = Progress(size(shock_range,1), dt = 1, desc="Creating transition probabilities for: ")
	for sh in 1:size(shock_range,1)
		row_range = (1+(sh-1)*endog_nstates*size(shock_range,1)):(sh)*endog_nstates*size(shock_range,1)
		#display(row_range)
		shocks_t = df_grid_new[row_range[1],exog_current]
		for vv in 1:len(exog_lead)
			transitions = transition_probabilities[exog_current[vv]]
			poss_trans = transitions[(transitions[:,1].==shocks_t[exog_current[vv]]),:]
			prob_densities[row_range,exog_lead_probs[vv]] = repeat(poss_trans[:,:prob],
				inner = prod(range_len[1:vv-1]), outer = endog_nstates*prod(range_len[(vv+1):nshocks]))
		end
		next!(prog)
	end

	df_grid_new.trans_prob = vec(prod(Matrix(prob_densities[:,exog_lead_probs]), dims = 2))::Array{Float64,1}

	# Check that the transition probs sum to the total number of initial points
	@assert round(sum(df_grid_new.trans_prob)) == round(nrow(df_grid)) "Transition probabilities don't sum to number of initial points"

	return df_grid_new
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
