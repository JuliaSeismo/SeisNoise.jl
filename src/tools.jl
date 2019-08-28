export snr, smooth, smooth!, nextpow2, abs_max!, standardize!, mad, savitsky_golay,
     running_mean, pws, std_threshold

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
    smooth(x,half_win=3)

Smooth Array `A` with half-window `half_win` (defaults to 3).
"""
function smooth!(A::AbstractArray; half_win::Int=3, dims::Int=1)
    # only use odd numbers
    if half_win % 2 != 1
        half_win += oneunit(half_win)
    end
    window_len = 2 * half_win + 1

    # get type of input Array
    T = eltype(A)

    if ndims(A) == 1
        A = reshape(A,size(A)...,1) # if 1D array, reshape to (length(A),1)
    end

    # extending the data at beginning and at the end
    # to apply the window at the borders
    X = vcat(A[window_len:-1:2,:],A,A[end-1:-1:end-window_len+1,:])

    # convolve with boxcar
    w = rect(window_len)
    w ./= sum(w)
    w = T.(w)
    for ii = 1:size(X,2)
        A[:,ii] .= conv(X[:,ii],w)[window_len+half_win:end-window_len-half_win+1]
    end
    return nothing
end
smooth(A::AbstractArray;half_win::Int=3, dims::Int=1) =
      (U = deepcopy(A);smooth!(U,half_win=half_win,dims=dims);return U)

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
    running_mean(x,N)

Returns array `x` smoothed by running mean of `N` points.
If N is even, reduces N by 1.
"""
function running_mean(x::AbstractArray,N::Int)
    if N % 2 == 0
        N -= 1
    end
    halfWindow = div(N,2)
    paddedx = vcat(x[1]*ones(halfWindow), x, x[end]*ones(halfWindow))
    y = conv(paddedx,ones(N) / N)
    return  y[2*halfWindow+1:end-2*halfWindow]
end

"""
    pws(A,fs,power,timegate)

Performs phase-weighted stack on array `A` of time series.

Follows methods of Schimmel and Paulssen, 1997.
If s(t) is time series data,
S(t) = s(t) + i*H(s(t)), where H(s(t)) is Hilbert transform of s(t)
S(t) = s(t) + i*H(s(t)) = A(t)*exp(i*phi(t)), where
A(t) is envelope of s(t) and phi(t) is phase of s(t)
Phase-weighted stack, g(t), is then:
g(t) = 1/N sum j = 1:N s_j(t) * | 1/N sum k = 1:N exp[i * phi_k(t)]|^v
where N is number of traces used, v is sharpness of phase-weighted stack

"""
function pws(A::AbstractArray, fs::Float64; power::Int=2, timegate::Float64=1.)
    M,N = size(A)
    analytic = A .+ im .* hilbert(A)
    phase = angle.(analytic)
    phase_stack = mean(exp.(im.*phase),dims=2)[:,1] ./ N # reduce array dimension
    phase_stack = abs.(phase_stack).^power

    # smoothing
    timegate_samples = convert(Int,timegate * fs)
    phase_stack = running_mean(phase_stack,timegate_samples)
    weighted = mean(A .* phase_stack,dims=2)
    return weighted
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
