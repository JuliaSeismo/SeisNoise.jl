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
function detrend!(X::AbstractArray{Float64,1})
    N = length(X)
    A = ones(N,2)
    A[:,1] = Array(1:N) ./ N
    coeff = lstsq(A,X)
    X[:] = X .- A *coeff
    return nothing
end
detrend(A::AbstractArray{Float64,1}) = (U = deepcopy(A);detrend!(U);return U)

function detrend!(X::AbstractArray{Float64,2})
    M,N = size(X)
    A = ones(M,2)
    A[:,1] = Array(1:N) ./ N
    for ii = 1:N
        coeff = lstsq(A,X[:,ii])
        X[:,ii] = X[:,ii] .- A *coeff
    end
    return X
end
detrend(A::AbstractArray{Float64,2}) = (U = deepcopy(A);detrend!(U);return U)

"""
    demean!(A)

Remove mean from columns of array `A`.
"""
function demean!(A::AbstractArray{Float64,1})
      μ = mean(A)
      for ii = 1:length(A)
        A[ii] -= μ
      end
  return A
end
demean(A::AbstractArray{Float64,1}) = (U = deepcopy(A);demean!(U);return U)
​
function demean!(A::AbstractArray{Float64,2})
      M,N = size(A)
      for ii = 1:N
        μ = mean(A[:,ii])
        for jj = 1:M
          A[jj,ii] -= μ
        end
      end
  return A
end
demean(A::AbstractArray{Float64,2}) = (U = deepcopy(A);demean!(U);return U)

end
