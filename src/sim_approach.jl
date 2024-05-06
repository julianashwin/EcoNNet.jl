"""
Functions that run the simulation approach for neural network learning
"""



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
            #predictions = cat(Array(s[tt-1, indices.outputindex_current]), Array(s[tt-1,indices.outputindex_lead]), dims = 1);
        else
            inputs::Array{Float64,1} = extract_inputs(s,tt,indices,options);
            predictions = predict!(inputs, beliefs);
            # For numerical stability at beginning of simulation, bound the expectations
            true_preds[tt,:] = predictions
            predictions = max.(predictions, -1000)
            predictions = min.(predictions, 1000)
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
            if tt%options.plotting_gap == 0
                if tt > (options.plotting_window+1)
                    plot(legend = :bottomleft, xguidefontsize=8)
                    for vv in options.plot_vars
                        (plot!(s[(tt - options.plotting_window):tt,vv], label = string(vv), legend = :bottomleft,xguidefontsize=8))
                    end
                    display(plot!(title = ""))
                else
                    plot(legend = :bottomleft, xguidefontsize=8)
                    for vv in options.plot_vars
                        (plot!(s[(1:options.plotting_window),vv], label = string(vv), legend = :bottomleft,xguidefontsize=8))
                    end
                    display(plot!(title = ""))
                end
            end
        end
    end

    return beliefs,s
end











"""
Function that simulates RLS learning and returns simulated data and updated beliefs
"""
function simulate_learning(sim_range::UnitRange{Int64}, s::DataFrame, beliefs::Dict, indices::EcoNNetIndices, options::EcoNNetOptions)

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
            predictions = max.(predictions, -100)
            predictions = min.(predictions, 100)
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
            beliefs::Dict = learn!(beliefs, s, tt, options, indices, loss)
            #global beliefs = learn!(beliefs, s, tt, options, indices, loss, ss_reminder = DNWR)
        end

        # Plot the evolution of variables
        if options.show_plots
            if tt%options.plotting_gap == 0 && tt > (options.plotting_window+1)
                display(join(["Simulation ", tt]))
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
