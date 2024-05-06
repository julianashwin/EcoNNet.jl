"""
DHM test 
"""
function DHM_test(s, ranges, nsamples; hvars = [], include_const = true)
    DHM_stats = zeros(length(ranges))
    for (ii, last_period) in enumerate(ranges)
        range = (last_period-nsamples+1):last_period
        errors = s.π[range] - s.Eπ_lead[range.-1]
        # Define h as function of state variables
        if include_const & length(hvars) > 0
            h = hcat(ones(nsamples), Matrix(s[range.-1, hvars]))
        elseif include_const & length(hvars) == 0
            h = hcat(ones(nsamples))
        else 
            h = hcat(Matrix(s[range.-1, hvars]))
        end
        M = sum(h.*errors, dims = 1)./nsamples
        # Define W matrix
        W = errors[1]*h[1,:]*h[1,:]'errors[1]./nsamples
        for tt in 2:nsamples
            W .+= errors[tt]*h[tt,:]*h[tt,:]'errors[tt]./nsamples
        end
        J = (nsamples*M*inv(W)*transpose(M))
        @assert size(J) == (1,1)
        DHM_stats[ii] = J[1]
    end
    nh = size(h,2)

    pvals = 1 .- cdf.([Chisq(nh)],DHM_stats)
    nsmall = round(sum(pvals .<= 0.05)/length(pvals), digits = 3)
    nlarge = round(sum(pvals .>= 0.95)/length(pvals), digits = 3)
    print("DHM test outcome is ", string(nsmall), " below and ", string(nlarge), " above.")

    return pvals
end

