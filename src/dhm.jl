"""
DHM test 
"""
function DHM_test(s, ranges, nsamples; hvars = [], include_const = true)
    DHM_stats = zeros(length(ranges))
    H = Matrix{Float64}(undef, 0, 0)
    for (ii, last_period) in enumerate(ranges)
        range = (last_period-nsamples+1):last_period
        errors = s.π[range] - s.Eπ_lead[range.-1]
        # Define h as function of state variables
        if include_const & length(hvars) > 0
            H = hcat(ones(nsamples), Matrix(s[range.-1, hvars]))
        elseif include_const & length(hvars) == 0
            H = hcat(ones(nsamples))
        else 
            H = hcat(Matrix(s[range.-1, hvars]))
        end
        M = sum(H.*errors, dims = 1)./nsamples
        # Define W matrix
        W = errors[1]*H[1,:]*H[1,:]'errors[1]./nsamples
        for tt in 2:nsamples
            W .+= errors[tt]*H[tt,:]*H[tt,:]'errors[tt]./nsamples
        end
        J = (nsamples*M*inv(W)*transpose(M))
        @assert size(J) == (1,1)
        DHM_stats[ii] = J[1]
    end
    nh = size(H,2)

    pvals = 1 .- cdf.([Chisq(nh)],DHM_stats)
    nsmall = round(sum(pvals .<= 0.05)/length(pvals), digits = 3)
    nlarge = round(sum(pvals .>= 0.95)/length(pvals), digits = 3)
    print("DHM test outcome is ", string(nsmall), " below and ", string(nlarge), " above.")

    return pvals
end

