export process_raw, process_raw!, process_fft, compute_fft, whiten, save_fft
import ArrayFuncs
"""
    compute_fft(S::SeisData,fs::Float64,freqmin::Float64,freqmax::Float64,
                         cc_step::Int, cc_len::Int;
                         time_norm::Union{Bool,String}=false,
                         to_whiten::Bool=false)


Computes windowed fft of ambient noise data.

Cross-correlates data using either cross-correlation, deconvolution or
cross-coherence. Saves cross-correlations in JLD2 data set.

# Arguments
- `S::SeisChannel`: SeisData structure.
- `fs::Float64`: Sampling rate to downsample `S`.
- `freqmin::Float64`: minimun frequency for whitening
- `freqmax::Float64`: maximum frequency for whitening
- `cc_step::Int`: time, in seconds, between successive cross-correlation windows
- `cc_len::Int`: length of noise data window, in seconds, to cross-correlate
- `time_norm::Union{Bool,String}`: time domain normalization to perform
- `to_whiten::Bool`: Apply whitening in frequency domain
"""
function compute_fft(S::SeisData,fs::Float64,freqmin::Float64,freqmax::Float64,
                     cc_step::Int, cc_len::Int;
                     time_norm::Union{Bool,String}=false,
                     to_whiten::Bool=false)

    # sync!(S,s=starttime,t=endtime)
    process_raw!(S,fs)  # demean, detrend, taper, lowpass, downsample
    merge!(S)
    starttime, endtime = u2d.(nearest_start_end(S[1],cc_len, cc_step))
    sync!(S,s=starttime,t=endtime)
    A, starts, ends = slide(S[1], cc_len, cc_step)
    FFT = process_fft(A, freqmin, freqmax, fs, time_norm=time_norm,
                      to_whiten=to_whiten)
    return F = FFTData(S[1].name, Dates.format(u2d(starts[1]),"Y-mm-dd"),
                       S[1].loc, S[1].fs, S[1].gain, freqmin, freqmax,
                       cc_len, cc_step, to_whiten, time_norm, S[1].resp,
                       S[1].misc, S[1].notes, starts, FFT)
end

"""

    data_checks!(S::SeisData)

Perform sanity checks on input raw data.
"""
function data_checks!(S::SeisData)
end

"""
    process_raw!(S,fs)

Pre-process raw seismic data.

Checks:
- sample rate is fs
- downsamples data
- checks for gaps in data
- phase-shifts data to begin at 00:00:00.0

# Arguments
- `S::SeisChannel`: SeisData structure.
- `fs::Float64`: Sampling rate to downsample `S`.
"""
function process_raw!(S::SeisData, fs::Float64)
    demean!(S)        # remove mean from channel
    ungap!(S)         # replace gaps with mean of channel
    if fs ∉ S.fs
        detrend!(S)       # remove linear trend from channel
        taper!(S)         # taper channel ends
        lowpass!(S,fs/2)    # lowpass filter before downsampling
        S = downsample(S,fs) # downsample to lower fs
    end
    phase_shift!(S) # timing offset from sampling period
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
- `freqmin::Float64`: minimun frequency for whitening
- `freqmax::Float64`: maximum frequency for whitening
- `time_norm::Union{Bool,String}`: time domain normalization to perform
- `to_whiten::Bool`: Apply whitening in frequency domain
- `corners::Int`: Number of corners in Butterworth filter.
- `zerophase::Bool`: If true, apply Butterworth filter twice for zero phase
                     change in output signal.
"""
function process_fft(A::AbstractArray,fs::Float64,freqmin::Float64,
                     freqmax::Float64;
                     time_norm::Union{Bool,String}=false,
                     to_whiten::Bool=false,
                     corners::Int=4,
                     zerophase::Bool=true)

    window_samples, N = size(A)

    # pre-process each window
    ArrayFuncs.demean!(A)
    ArrayFuncs.detrend!(A)
    ArrayFuncs.taper!(A,fs)
    ArrayFuncs.bandpass!(A,freqmin,freqmax,fs,corners=corners,
                         zerophase=zerophase)
    ArrayFuncs.demean!(A)

    # apply one-bit normalization
    if time_norm == "one_bit"
        A .= sign.(A)
    end

    # take fft along first dimension
    if to_whiten
        FFT = whiten(A,freqmin, freqmax, fs)
    else
        FFT = rfft(A,1)
    end
    return FFT
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
    remove_resp(args)

remove instrument response - will require reading stationXML and extracting poles
and zeros
"""
function remove_resp(args)
end

"""
    save_fft(F::FFTData, OUT::String)

Save FFTData `F` to JLD2.
"""
function save_fft(F::FFTData, FFTOUT::String)
    # check if FFT DIR exists
    if isdir(FFTOUT) == false
        mkpath(FFTOUT)
    end

    # create JLD2 file and save FFT
    net,sta,loc,chan = split(F.name,'.')
    filename = joinpath(FFTOUT,"$net.$sta.jld2")
    file = jldopen(filename, "a+")
    if !(chan in keys(file))
        group = JLD2.Group(file, chan)
        group[F.id] = F
    else
        file[chan][F.id] = F
    end
    close(file)
end
