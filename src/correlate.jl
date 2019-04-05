module Correlate
# cross-correlation module
export clean_up, lstsq, detrend, taper, whiten, correlate
using ..SeisJul

"""
    clean_up(A,fs,freqmin,freqmax)

Demean, detrend, taper and filter time series.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
"""
function clean_up(A::AbstractArray, fs::Real, freqmin::Real, freqmax::Real)
    clean = A .- mean(A,dims=1)
    clean = detrend(clean)
    clean = taper(clean,fs)
    clean = bandpass(clean,freqmin,freqmax,fs)
    return clean
end

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
    detrend(A)

Remove linear trend from array `A` using least-squares regression.
"""
function detrend(X::AbstractArray)
    N = length(X)
    A = ones(N,2)
    A[:,1] = Array(1:N) ./ N
    coeff = lstsq(A,X)
    newX = X .- A *coeff
end

"""
    taper(A,fs; max_percentage=0.05, type="hann", max_length=20)

Taper a time series `A` with sampling_rate `fs`.
Defaults to 'hann' window. Uses smallest of `max_percentage` * `fs`
or `max_length`.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
- `max_percentage::float`: Decimal percentage of taper at one end (ranging
    from 0. to 0.5).
- `type::String`: Type of taper to use for detrending. Defaults to 'hann'.
- `max_length::Real`: Length of taper at one end in seconds.
"""
function taper(A::AbstractArray, fs::Real; max_percentage=0.05,
               type="hann", max_length::Real=20)
    N = length(A)
    wlen = min(Int(N * max_percentage), Int(max_length * fs), Int(N/2))
    taper_sides = -hanning(2 * wlen, zerophase=true) .+ 1
    tape = vcat(taper_sides[1:wlen], ones(N - 2 * wlen), taper_sides[wlen+1:end])
    newA = A .* tape
end

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
function whiten(A::AbstractArray, fs::Real, freqmin::Real, freqmax::Real; pad::Int=100)
    N = length(A)

    # get whitening frequencies
    freqvec = rfftfreq(N,fs)
    freqind = findall(x -> x >= freqmin && x <= freqmax, freqvec)
    low, high = freqind[1] - pad, freqind[end] + pad
    left, right = freqind[1], freqind[end]

    if low <= 1
        low = 1
    end

    if high > length(freqvec) / 2
        high = Int(length(freqind) / 2)
    end

    # take fft and whiten
    fftraw = rfft(A,1)
    # left zero cut-off
    fftraw[1:low] .= 0. + 0.0im
    # left tapering
    fftraw[low+1:left] .= cos.(LinRange(pi / 2., pi, left - low)).^2 .* exp.(
        im .* angle.(fftraw[low:left-1]))
    # pass band
    fftraw[left:right] .= exp.(im .* angle.(fftraw[left:right]))
    # right tapering
    fftraw[right+1:high] .= cos.(LinRange(0., pi/2., high-right)).^2 .* exp.(
        im .* angle.(fftraw[right+1:high]))
    # right zero cut-off
    fftraw[high+1:end] .= 0. + 0.0im
    return fftraw
end

"""
    correlate(fft1, fft2, maxlag, method='cross-correlate')

Cross-correlation of two ffts.


"""
function correlate(fft1::AbstractArray, fft2::AbstractArray, N::Int, maxlag::Int;
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

end
