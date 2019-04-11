# cross-correlation module
export clean_up!, clean_up, whiten, correlate, next_fast_len

"""
    clean_up!(A,fs,freqmin,freqmax)

Demean, detrend, taper and filter time series.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
"""
function clean_up!(A::AbstractArray, fs::Real, freqmin::Real, freqmax::Real;
                   corners::Int=4, zerophase::Bool=false)
    ArrayFuncs.demean!(A)
    ArrayFuncs.detrend!(A)
    taper!(A,fs)
    bandpass!(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
    return nothing
end
clean_up(A::AbstractArray, fs::Real, freqmin::Real, freqmax::Real) =
                 (U = deepcopy(A); clean_up!(U,fs,freqmin,freqmax); return U)



"""
    whiten(A, fs, freqmin, freqmax, pad=100)

Whiten spectrum of time series `A` between frequencies `freqmin` and `freqmax`.
Uses real fft to speed up computation.
Returns the whitened (single-sided) fft of the time series.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
- `pad::Int`: Number of tapering points outside whitening band.
"""
function whiten(A::AbstractArray, freqmin::Real, freqmax::Real, fs::Real;
                pad::Int=100)
    if ndims(A) == 1
        A = reshape(A,size(A)...,1) # if 1D array, reshape to (length(A),1)
    end

    N,_ = size(A)

    # get whitening frequencies
    freqvec = rfftfreq(N,fs)
    freqind = findall(x -> x >= freqmin && x <= freqmax, freqvec)
    low, high = freqind[1] - pad, freqind[end] + pad
    left, right = freqind[1], freqind[end]

    if low <= 1
        low = 1
        left = low + 100
    end

    if high > length(freqvec)
        high = length(freqvec)- 1
        right = high - 100
    end

    # take fft and whiten
    fftraw = rfft(A,1)
    # left zero cut-off
    fftraw[1:low,:] .= 0. + 0.0im
    # left tapering
    fftraw[low+1:left,:] .= cos.(LinRange(pi / 2., pi, left - low)).^2 .* exp.(
        im .* angle.(fftraw[low:left-1,:]))
    # pass band
    fftraw[left:right,:] .= exp.(im .* angle.(fftraw[left:right,:]))
    # right tapering
    fftraw[right+1:high,:] .= cos.(LinRange(0., pi/2., high-right)).^2 .* exp.(
        im .* angle.(fftraw[right+1:high,:]))
    # right zero cut-off
    fftraw[high+1:end,:] .= 0. + 0.0im
    return fftraw
end

"""
    correlate(fft1, fft2, maxlag, method='cross-correlate')

Cross-correlation of two ffts.


"""
function cross_corr(fft1::AbstractArray, fft2::AbstractArray, N::Int, maxlag::Int;
                   method::String="cross-correlation")

    corrF = conj.(fft1) .* fft2
    if method == "deconv"
        corrF ./= abs.(fft1).^2
    elseif method == "coherence"
        corrF ./= abs.(fft1)
        corrF ./= abs.(fft2)
    end

    # take inverse fft
    corrT = irfft(corrF,N)
    corrT = fftshift(corrT)

    # return corr[-maxlag:maxlag]
    t = range(-Int(N/2) + 1, Int(N/2) - 1)
    ind = findall(x -> abs(x) <= maxlag,t)
    corrT = corrT[ind]
end

"""
    next_fast_len(N::Real)


"""
function next_fast_len(N::Real)
    return nextprod([2,3,5],N)
end
