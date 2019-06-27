export process_raw, process_raw!, process_fft, compute_fft, whiten
"""
    compute_fft(S::SeisData,freqmin::Float64,freqmax::Float64, fs::Float64
                         cc_step::Int, cc_len::Int;
                         time_norm::Union{Bool,String}=false,
                         to_whiten::Bool=false)


Computes windowed fft of ambient noise data.

Cross-correlates data using either cross-correlation, deconvolution or
cross-coherence. Saves cross-correlations in JLD2 data set.

# Arguments
- `S::SeisChannel`: SeisData structure.
- `freqmin::Float64`: minimun frequency for whitening.
- `freqmax::Float64`: maximum frequency for whitening.
- `fs::Float64`: Sampling rate to downsample `S`.
- `cc_step::Int`: time, in seconds, between successive cross-correlation windows.
- `cc_len::Int`: length of noise data window, in seconds, to cross-correlate.
- `time_norm::Union{Bool,String}`: time domain normalization to perform.
- `to_whiten::Bool`: Apply whitening in frequency domain.
"""
function compute_fft(S::SeisData,freqmin::Float64,freqmax::Float64,fs::Float64,
                     cc_step::Int, cc_len::Int;
                     time_norm::Union{Bool,String}=false,
                     to_whiten::Bool=false)

    # sync!(S,s=starttime,t=endtime)
    merge!(S)
    ungap!(S)
    process_raw!(S,fs)  # demean, detrend, taper, lowpass, downsample
    starttime, endtime = u2d.(nearest_start_end(S[1],cc_len, cc_step))
    sync!(S,s=starttime,t=endtime)
    A, starts, ends = slide(S[1], cc_len, cc_step)
    FFT = process_fft(A, freqmin, freqmax, fs, time_norm=time_norm,
                      to_whiten=to_whiten)
    return F = FFTData(S[1].id, Dates.format(u2d(starts[1]),"Y-mm-dd"),
                       S[1].loc, S[1].fs, S[1].gain, freqmin, freqmax,
                       cc_len, cc_step, to_whiten, time_norm, S[1].resp,
                       S[1].misc, S[1].notes, starts, FFT)
end

"""

    data_checks!(S::SeisData)

Perform sanity checks on input raw data.
"""
function data_checks!(S::SeisData)
    for ii = 1:S.n
        if S[ii].fs < 1
            S[ii].fs = 1 / S[ii].fs
        end
    end
end

"""
    process_raw!(S,fs)

Pre-process raw seismic data.

- Removes mean from each channel in `S`.
- Detrends each channel in `S`.
- Downsamples data to sampling rate `fs`
- Phase-shifts data to begin at 00:00:00.0

# Arguments
- `S::SeisData`: SeisData structure.
- `fs::Float64`: Sampling rate to downsample `S`.
"""
function process_raw!(S::SeisData, fs::Float64)
    for ii = 1:S.n
        demean!(S[ii].x)        # remove mean from channel
        if fs ∉ S.fs
            detrend!(S[ii].x)       # remove linear trend from channel
            taper!(S[ii].x,S[ii].fs)         # taper channel ends
            lowpass!(S[ii].x,fs/2,S[ii].fs)    # lowpass filter before downsampling
            S[ii] = downsample(S[ii],fs) # downsample to lower fs
        end
        phase_shift!(S[ii]) # timing offset from sampling period
    end
    return nothing
end
process_raw(S::SeisData, fs::Float64) = (U = deepcopy(S);
            process_raw!(U,fs); return U)

"""
    process_fft(A::AbstractArray,freqmin::Float64,freqmax::Float64,fs::Float64;
                time_norm=false,to_whiten=false,corners=corners,
                zerophase=zerophase)

# Arguments
- `A::AbstractArray`: Array with time domain data.
- `fs::Float64`: Sampling rate of data in `A`.
- `freqmin::Float64`: minimum frequency for whitening.
- `freqmax::Float64`: maximum frequency for whitening.
- `time_norm::Union{Bool,String}`: time domain normalization to perform.
- `to_whiten::Bool`: Apply whitening in frequency domain.
- `corners::Int`: Number of corners in Butterworth filter.
- `zerophase::Bool`: If true, apply Butterworth filter twice for zero phase
                     change in output signal.
"""
function process_fft(A::AbstractArray,freqmin::Float64, freqmax::Float64,
                     fs::Float64; time_norm::Union{Bool,String}=false,
                     to_whiten::Bool=false,
                     corners::Int=4,
                     zerophase::Bool=true)

    # pre-process each window
    demean!(A)
    detrend!(A)
    taper!(A,fs)
    bandpass!(A,freqmin,freqmax,fs,corners=corners,
                         zerophase=zerophase)
    demean!(A)

    if eltype(A) != Float32
        A = Float32.(A)
    end

    if to_whiten && time_norm == false
        M,N = size(A)
        FFT = whiten(A,freqmin, freqmax, fs)
    elseif to_whiten && time_norm != false
        M,N = size(A)
        FFT = whiten(A,freqmin, freqmax, fs)
        A = irfft(FFT,M,1)
        # apply time-domain normalization or extract instantaneous phase
        if time_norm == "one_bit"
            A .= sign.(A)
        elseif time_norm == "phase"
            phase!(A)
        end
        FFT = rfft(A,1)
    else
        if time_norm == "one_bit"
            A .= sign.(A)
        elseif time_norm == "phase"
            phase!(A)
        end
        FFT = rfft(A,1)
    end
    return FFT
end

"""
    whiten(A, freqmin, freqmax, fs, pad=100)

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
function whiten(A::AbstractArray, freqmin::Float64, freqmax::Float64, fs::Float64;
                pad::Int=50)
    if ndims(A) == 1
        A = reshape(A,size(A)...,1) # if 1D array, reshape to (length(A),1)
    end

    # get size and convert to Float32
    M,N = size(A)
    if eltype(A) != Float32
        A = Float32.(A)
    end

    # get whitening frequencies
    freqvec = rfftfreq(M,fs)
    left = findfirst(x -> x >= freqmin, freqvec)
    right = findlast(x -> x <= freqmax, freqvec)
    low, high = left - pad, right + pad

    if low <= 1
        low = 1
        left = low + pad
    end

    if high > length(freqvec)
        high = length(freqvec)- 1
        right = high - pad
    end

    # take fft and whiten
    fftraw = rfft(A,1)
    # left zero cut-off
     for jj = 1:N
         for ii = 1:low
            fftraw[ii,jj] = 0. + 0.0im
        end
    end
    # left tapering
     for jj = 1:N
         for ii = low+1:left
            fftraw[ii,jj] = cos(pi / 2 + pi / 2* (ii-low-1) / pad).^2 * exp(
            im * angle(fftraw[ii,jj]))
        end
    end
    # pass band
     for jj = 1:N
         for ii = left:right
            fftraw[ii,jj] = exp(im * angle(fftraw[ii,jj]))
        end
    end
    # right tapering
     for jj = 1:N
         for ii = right+1:high
            fftraw[ii,jj] = cos(pi/2 * (ii-right) / pad).^2 * exp(
            im * angle(fftraw[ii,jj]))
        end
    end
    # right zero cut-off
     for jj = 1:N
         for ii = high+1:size(fftraw,1)
            fftraw[ii,jj] = 0. + 0.0im
        end
    end
    return fftraw
end
