"""
Functions to visualise results
"""


"""
Turn number into string, replacing decimal points with the letter "p"
	This is useful for filenames that give some details on parameter values
"""
function rep_pnt(num)
	num_out = replace(string(num), "."=> "p")
	return num_out
end


"""
Function that generates a perfect foresight path from given starting point
"""
function pf_path(initial_ss; periods = 100, lags::Int64 = 1)

	paths = DataFrame(zeros(periods, len(variables)), variables);
    paths = initialise_df(paths, initial_ss);

    for tt in (lags+1):periods
        #print(tt)
		if lags == 1
			inputs = vcat(Array(paths[tt-1,indices.endognames]))
		elseif lags == 2
        	inputs = vcat(Array(paths[tt-2,indices.endognames]),Array(paths[tt-1,indices.endognames]))
		else
			display("More than 2 lags not supported")
		end
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
				plot!(paths[:,plot_vars[vv]], label = string(plot_vars[vv]))
			else
				display(plot!(paths[:,plot_vars[vv]], label = string(plot_vars[vv])))
			end
		end
	end

    return paths

end


function irf(type, initial_ss, beliefs::Dict; periods = 100, magnitude = 1.0, persistence = 0.9,
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
	h_points = 11:99, v_points = 12:100, label = "", arrow_size = .4,
	final_arrow = false)

	h_path = paths[h_points,vars[1]]
	v_path =paths[v_points,vars[2]]

	if final_arrow
		plt = plot!(h_path,v_path,color = :blue, arrow = (arrow_size,arrow_size),label = label)
	else
		plt = plot!(h_path,v_path,color = :blue, label = label)
	end
	for point in arrow_points
		plot!([h_path[point],h_path[point+1]], [v_path[point],v_path[point+1]],
			arrow = (arrow_size,arrow_size), color = :blue, label = "")
	end
	display(plot!())
	return plt
end
