export snr, smooth, smooth!, nextpow2, abs_max!, standardize!, mad, savitsky_golay,
     movingaverage, pws, std_threshold

"""
    snr(A,fs)

Signal to noise ratio of cross-correlations in matrix `A` with `fs'.

Follows method of Clarke et. al, 2011. Measures SNR at each point.
"""
function snr(A::AbstractArray,fs::Real)
    t,N = size(A)
    A_mean = mean(A,dims=2)

    # calculate noise and envelope functions
    sigma = mean(A.^2,dims=2) .- A_mean.^2
    sigma = sqrt.(sigma./ (N-1))
    return sigma
end
snr(C::CorrData) = snr(C.corr,C.fs)

"""
    smooth(A, half_win)

Smooth array `A` with half-window `half_win` (defaults to 3).
"""
function smooth!(A::AbstractArray, half_win::Int=3)
    if ndims(A) == 1
        return movingaverage(A,half_win)
    end

    Nrows, Ncols = size(A)

    for ii = 1:Ncols
        A[:,ii] .= movingaverage(A[:,ii],half_win)
    end
    return nothing
end
smooth(A::AbstractArray,half_win::Int=3) =
      (U = deepcopy(A);smooth!(U,half_win);return U)

"""
  movingaverage(A,half_win)

Smooth array `A` with moving average of half-window `half-win` (defaults to 3.)
"""
function movingaverage(A::AbstractArray, half_win::Int=3)
    prepend!(A,A[1:half_win])
    append!(A,A[end-half_win:end])
    N = length(A)
    B = zeros(eltype(A),N)
    window_len = 2 * half_win + 1
    s = sum(A[1:window_len])
    B[half_win+1] = s
    for ii = half_win+2:N-half_win
        s = s - A[ii-half_win] + A[ii+half_win]
        B[ii] = s
    end
    B ./= window_len
    return B[half_win+1:end-half_win-1]
end

"""
    nextpow2(x)

Returns the next power of 2 of real, positive number `x`.
"""
function nextpow2(x::Real)
    ceil(Int,log2(abs(x)))
end

"""
    abs_max!(A)

Returns array `A` divided by its absolute maximum value.
"""
function abs_max!(A::AbstractArray)
    A .= A ./ maximum(abs.(A),dims=1)
    return nothing
end
abs_max!(C::CorrData) = abs_max!(C.corr)

"""
    standardize!(A)

Demean and standardize array `A` to unit std.
"""
function standardize!(A::AbstractArray)
    A .= (A .- mean(A,dims=1)) ./ std(A,dims=1)
    return nothing
end
standardize!(C::CorrData) = standardize!(C.corr)

"""
    mad(A)

Median Absolute Deviation of array `A`.
MAD = median(|Xi- median(A)|)
"""
function mad(A::AbstractArray)
    return median(abs.(A .- median(A,dims=1)),dims=1)
end
mad(C::CorrData) = mad(C.corr)

"""
    savitsky_golay(x, window, polyOrder, [deriv])

Polynomial smoothing of vector `x` with Savitsky Golay filter.
Polynomial order `polyOrder' must be less than window length `N`.
Modified from https://github.com/BBN-Q/Qlab.jl/blob/master/src/SavitskyGolay.jl
"""
function savitsky_golay(x::Vector, N::Int, polyOrder::Int; deriv::Int=0)
    #Some error checking
    @assert isodd(N) "Window size must be an odd integer."
    @assert polyOrder < N "Polynomial order must me less than window size."

    halfWindow = Int((N-1)/2)

    #Setup the S matrix of basis vectors.
    S = zeros(N, polyOrder+1)
    for ct = 0:polyOrder
        S[:,ct+1] = Array(-halfWindow:halfWindow).^ct
    end

    #Compute the filter coefficients for all orders
    G = S*pinv(S'*S)

    #Slice out the derivative order we want
    filterCoeffs = G[:,deriv+1] * factorial(deriv);

    #Pad the signal with the endpoints and convolve with filter
    paddedX = vcat(x[1]*ones(halfWindow), x, x[end]*ones(halfWindow))
    y = conv(filterCoeffs[end:-1:1], paddedX)

    #Return the valid midsection
    return y[2*halfWindow+1:end-2*halfWindow]
end

"""

  std_threshold(A,thresh)

Returns indices of cols of `A` where max(abs.A[:,col]))/std(A[:,col]) < `max_std`.
"""
function std_threshold(A::AbstractArray, max_std::Float64)
    stds = std(A,dims=1)[1,:]
    maxs = maximum(abs.(A),dims=1)[1,:]
    threshs = maxs ./ stds
    ind = []
    for ii = 1:length(threshs)
        if threshs[ii] < max_std
            append!(ind,ii)
        end
    end
    return ind
end
