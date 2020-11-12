"""
Functions for the state variable dataframe used
"""


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
