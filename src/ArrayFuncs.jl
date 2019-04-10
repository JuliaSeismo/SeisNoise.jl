module ArrayFuncs
import Statistics.mean
import LinearAlgebra.pinv
"""
    lstsq(A,X)

Least-squares regression of array `A` using the pseudo-inverse.

Solves the equation `A X = B` by computing a vector `X` that
    minimizes the Euclidean 2-norm `|| B - A X ||^2`.

# Arguments
- `A::AbstractArray`: Coefficient matrix.
- `X::AbstractArray`: Dependent variable.
"""
function lstsq(A::AbstractArray,X::AbstractArray)
    coeff = pinv(A' * A) * A' * X
end

"""
    detrend!(A)

Remove linear trend from array `A` using least-squares regression.
"""
function detrend!(X::AbstractArray)
    N = length(X)
    A = ones(N,2)
    A[:,1] = Array(1:N) ./ N
    coeff = lstsq(A,X)
    X[:] = X .- A *coeff
    return nothing
end
detrend(A::AbstractArray) = (U = deepcopy(A);detrend!(U);return U)

"""
    demean!(A)

Remove mean from array `A`.
"""
function demean!(A::AbstractArray)
    A[:] = A .- mean(A,dims=1)
    return nothing
end
demean(A::AbstractArray) = (U = deepcopy(A);demean!(U);return U)

end
