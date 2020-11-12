"""
Auxilliary functions to make code more elegant etc...
"""

"""
A shortened name for the length function to make code more compact
"""
function len(x)
	length(x)
end



"""
A function that solves a quadratic equation
"""
function solvequadratic(a, b, c)
    d = sqrt(Complex(b^2 - 4a*c))
    solns = (-b - d) / 2a, (-b + d) / 2a
    if d != real.(d)
        println("No real solutions")
    else
        solns = real.(solns)
    end
	if a == 0.0
		solns = -c/b
	elseif c == 0.0
		solns = -b/a
	end
    return solns

end


"""
Function that finds the nearest value to x in an array
	A is an array of options from which the nearest will be returned
	x is the array which needs to be matched to a row of A
"""
function findnearest(A::Array{Float64,1},x::Float64)

	indx = findmin(abs.(A.-x))[2]
	soln = A[indx]

	return soln
end
