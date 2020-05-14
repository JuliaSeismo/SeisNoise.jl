export snr, smooth, smooth!, abs_max, abs_max!, standardize, standardize!
export mad, std_threshold

"""
    snr(A)

Signal to noise ratio of cross-correlations in matrix `A`.

Follows method of Clarke et. al, 2011. Measures SNR at each point.
"""
function snr(A::AbstractArray)
    Nrows,Ncols = size(A)
    A_mean = mean(A,dims=2)

    # calculate noise and envelope functions
    sigma = mean(A.^2,dims=2) .- A_mean.^2
    sigma .= sqrt.(sigma./ (Ncols-1))
    s = abs.(A_mean .+ im .* hilbert(A_mean))
    return s ./ sigma
end
snr(C::CorrData) = snr(C.corr)

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
    abs_max!(A)

Returns array `A` divided by its absolute maximum value.
"""
function abs_max!(A::AbstractArray)
    maxabs = maximum(abs.(A),dims=1)
    if any(maxabs .== 0)
        throw(DomainError("All zero column leads to NaN"))
    end
    A ./= maxabs
    return nothing
end
abs_max(A::AbstractArray) = (U = deepcopy(A);abs_max!(U);return U)
abs_max!(C::CorrData) = abs_max!(C.corr)
abs_max(C::CorrData) = (U = deepcopy(C);abs_max!(U);return U)

"""
    standardize!(A)

Demean and standardize array `A` to unit std.
"""
function standardize!(A::AbstractArray)
    A .-= mean(A,dims=1)
    A ./= std(A,dims=1)
    return nothing
end
standardize(A::AbstractArray) = (U = deepcopy(A);standardize!(U);return U)
standardize!(C::CorrData) = standardize!(C.corr)
standardize(C::CorrData) = (U = deepcopy(C);standardize!(U);return U)

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
