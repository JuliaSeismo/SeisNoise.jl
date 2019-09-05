export process_raw, process_raw!, process_fft, compute_fft, whiten, remove_response!, read_stationXML
"""
    compute_fft(S, freqmin, freqmax, fs, cc_step, cc_len;
                time_norm=false, to_whiten=false, max_std=5.)


Computes windowed fft of ambient noise data.

# Arguments
- `S::SeisData`: SeisData structure.
- `freqmin::Float64`: minimum frequency for filtering/whitening.
- `freqmax::Float64`: maximum frequency for filtering/whitening.
- `fs::Float64`: Sampling rate to downsample `S`.
- `cc_step::Int`: time, in seconds, between successive cross-correlation windows.
- `cc_len::Int`: length of noise data window, in seconds, to cross-correlate.
- `time_norm::Union{Bool,String}`: time domain normalization to perform.
- `to_whiten::Bool`: Apply whitening in frequency domain.
- `max_std::Float64=5.`: Number of standard deviations above mean to reject windowed data.
"""
function compute_fft(S::SeisData,freqmin::Float64,freqmax::Float64,fs::Float64,
                     cc_step::Int, cc_len::Int;
                     time_norm::Union{Bool,String}=false,
                     to_whiten::Bool=false,
                     max_std::Float64=5.)

    # sync!(S,s=starttime,t=endtime)
    merge!(S)
    ungap!(S)
    process_raw!(S,fs)  # demean, detrend, taper, lowpass, downsample
    starttime, endtime = u2d.(nearest_start_end(S[1],cc_len, cc_step))
    sync!(S,s=starttime,t=endtime)
    A, starts, ends = slide(S[1], cc_len, cc_step)
    ind = std_threshold(A,max_std)
    if length(ind) == 0
        error("No windows remaining for day $(Dates.format(u2d(starts[1]),"Y-mm-dd"))
              with max std = $max_std.")
    end
    A = A[:,ind]
    starts = starts[ind]
    FFT = process_fft(A, freqmin, freqmax, fs, time_norm=time_norm,
                      to_whiten=to_whiten)
    return F = FFTData(S[1].id, Dates.format(u2d(starts[1]),"Y-mm-dd"),
                       S[1].loc, S[1].fs, S[1].gain, freqmin, freqmax,
                       cc_len, cc_step, to_whiten, time_norm, S[1].resp,
                       S[1].misc, S[1].notes, starts, FFT)
end

"""
    compute_fft(S, freqmin, freqmax, fs, cc_step, cc_len, stationXML,
                time_norm=false,to_whiten=false, max_std=5.)


Computes windowed fft of ambient noise data.

Removes instrument response from `S` using response from `stationXML`.

# Arguments
- `S::SeisData`: SeisData structure.
- `freqmin::Float64`: minimum frequency for instrument response pre-filter.
- `freqmax::Float64`: maximum frequency for instrument response pre-filter.
- `fs::Float64`: Sampling rate to downsample `S`.
- `cc_step::Int`: time, in seconds, between successive cross-correlation windows.
- `cc_len::Int`: length of noise data window, in seconds, to cross-correlate.
- `time_norm::Union{Bool,String}=false`: time domain normalization to perform.
- `to_whiten::Bool=false`: Apply whitening in frequency domain.
- `max_std::Float64`: Number of standard deviations above mean to reject windowed data.
"""
function compute_fft(S::SeisData,freqmin::Float64,freqmax::Float64,fs::Float64,
                     cc_step::Int, cc_len::Int, stationXML::String;
                     time_norm::Union{Bool,String}=false,
                     to_whiten::Bool=false, max_std::Float64=5.,
                     whitemin::Union{Nothing,Float64}=nothing,
                     whitemax::Union{Nothing,Float64}=nothing)

    # sync!(S,s=starttime,t=endtime)
    merge!(S)
    ungap!(S)
    process_raw!(S,fs)  # demean, detrend, taper, lowpass, downsample
    starttime, endtime = u2d.(nearest_start_end(S[1],cc_len, cc_step))
    sync!(S,s=starttime,t=endtime) # sync start and end times
    remove_response!(S,stationXML,freqmin,freqmax) # remove inst response
    A, starts, ends = slide(S[1], cc_len, cc_step) # cut waveform into windows

    # amplitude threshold indices
    ind = std_threshold(A,max_std)
    if length(ind) == 0
        error("No windows remaining for day $(Dates.format(u2d(starts[1]),"Y-mm-dd"))
              with max std = $max_std.")
    end
    A = A[:,ind]
    starts = starts[ind]

    # check whitening frequencies
    if isnothing(whitemin)
        whitemin = freqmin
    end
    if isnothing(whitemax)
        whitemax = freqmax
    end
    FFT = process_fft(A, whitemin, whitemax, fs, time_norm=time_norm,
                      to_whiten=to_whiten)
    return F = FFTData(S[1].id, Dates.format(u2d(starts[1]),"Y-mm-dd"),
                       S[1].loc, S[1].fs, S[1].gain, freqmin, freqmax,
                       cc_len, cc_step, to_whiten, time_norm, S[1].resp,
                       S[1].misc, S[1].notes, starts, FFT)
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

"""
  remove_response!(S, stationXML, freqmin, freqmax)

Loads instrument response from stationXML and removes response from `S`.

# Arguments
- `S::SeisData`: SeisData structure.
- `stationXML::String`: Path to stationXML file, e.g. "/path/to/file.xml"
- `freqmin::Float64`: minimum frequency for pre-filtering.
- `freqmax::Float64`: maximum frequency for pre-filtering.
- `np::Int`: number of poles in pre-filter.
- `t_max::Float64`: Length of taper in seconds.
- `wl::Float32`: Water level for instrument response.
"""
function remove_response!(S::SeisData, stationXML::String, freqmin::Float64,
                          freqmax::Float64;np::Int=2, t_max::Float64=20.,
                          wl::Float32=eps(Float32))

    # read response file
    if !isfile(stationXML)
        error("$stationXML does not exist. Instrument response not removed.")
    end

    R = read_stationXML(stationXML)
    # loop through responses
    Rid = R.id
    Sid = S.id

    # remove trend, taper and prefilter
    SeisIO.detrend!(S)
    SeisIO.taper!(S,t_max=t_max)
    filtfilt!(S,fl=freqmin,fh=freqmax,np=np,rt="Bandpass")

    # remove instrument response for each channel in S
    @inbounds for ii = 1:S.n
        id = S[ii].id
        ind = findfirst(x -> x == id,Rid)
        LOC = R[ind].loc
        GAIN = R[ind].gain
        RESP = R[ind].resp
        UNITS = R[ind].units
        setindex!(S.loc,LOC,ii)
        setindex!(S.gain,GAIN,ii)
        setindex!(S.units,UNITS,ii)
        translate_resp!(S,RESP,chans=ii,wl=wl)
    end
    return nothing
end
remove_response(S::SeisData, stationXML::String, freqmin::Float64,
                          freqmax::Float64;np::Int=2, t_max::Float64=20.,
                          wl::Float32=eps(Float32)) = (U = deepcopy(S);
            remove_response!(U,stationXML,freqmin,freqmax,np=np,t_max=t_max,
                             wl=wl); return U)

"""
  read_stationXML(stationXML)

Reads instrument response from stationXML file.

Returns a SeisData object with an instrument response for each channel in
the stationXML file.

# Arguments
- `stationXML::String`: Path to stationXML file, e.g. "/path/to/file.xml"

"""
function read_stationXML(stationXML::String)
    io = open(stationXML, "r")
    xsta = read(io, String)
    close(io)
    return SeisIO.FDSN_sta_xml(xsta)
end
