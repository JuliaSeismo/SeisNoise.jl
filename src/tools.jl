export snr, smooth, smooth!, nextpow2, abs_max!, standardize!, mad, savitsky_golay,
    pws, std_threshold

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
function smooth!(A::AbstractArray, half_win::Int=3, dims::Int=1)
    T = eltype(A)
    window_len = 2 * half_win + 1
    csumsize = tuple(collect(size(A)) .+ [i==1 for i in 1:ndims(A)]...)
    csum = similar(A,T,csumsize)
    csum[1,:] .= zero(T)
    csum[2:end,:] = cumsum(A,dims=dims)
    A[half_win+1:end-half_win,:] .= (csum[window_len+1:end,:] .- csum[1:end-window_len,:]) ./ window_len
    return nothing
end
smooth(A::AbstractArray,half_win::Int=3, dims::Int=1) =
      (U = deepcopy(A);smooth!(U,half_win,dims);return U)

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
