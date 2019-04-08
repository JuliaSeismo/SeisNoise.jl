export bandpass, bandpass!, bandstop, bandstop!, lowpass, lowpass!,
        highpass, highpass!, taper, taper!, envelope

"""
    bandpass!(A,freqmin,freqmax,fs,corners=4,zerophase=false)

Butterworth-Bandpass Filter.

Filter data `A` from `freqmin` to `freqmax` using `corners` corners.

# Arguments
- `A::AbstractArray`: Data to filter
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
- `fs::Real`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function bandpass!(A::AbstractArray, freqmin::Real, freqmax::Real, fs::Real;
                   corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    low = freqmin / fe
    high = freqmax / fe

    # warn if above Nyquist frequency
    if high - oneunit(high) > -1e-6
        @warn "Selected high corner frequency ($freqmax) of bandpass is at or
        above Nyquist ($fe). Applying a high-pass instead."
        highpass!(A,freqmin,fs,corners=corners,zerophase=zerophase)
        return nothing
    end

    # throw error if low above Nyquist frequency
    if low > 1
        ArgumentError("Selected low corner frequency is above Nyquist.")
    end

    # create filter
    responsetype = Bandpass(freqmin, freqmax; fs=fs)
    designmethod = Butterworth(corners)
    if zerophase
        A[:] = filtfilt(digitalfilter(responsetype, designmethod), A)
    else
        A[:] = filt(digitalfilter(responsetype, designmethod), A)
    end

    return nothing
end
bandpass(A::AbstractArray,freqmin::Real, freqmax::Real; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(A);bandpass!(A,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

function bandpass!(C::SeisChannel, freqmin::Real, freqmax::Real;
                   corners::Int=4, zerophase::Bool=false)
    fs = C.fs
    bandpass!(C.x,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
    return nothing
end
bandpass(C::SeisChannel,freqmin::Real, freqmax::Real; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(C);bandpass!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

function bandpass!(S::SeisData, freqmin::Real, freqmax::Real;
                   corners::Int=4, zerophase::Bool=false)
    @inbounds for i = 1:S.n
        bandpass!(S[i],freqmin,freqmax,corners=corners,zerophase=zerophase)
    end
    return nothing
end

bandpass(S::SeisData,freqmin::Real, freqmax::Real; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(S);bandpass!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

"""
    bandstop!(A,freqmin,freqmax,fs,corners=4,zerophase=false)

Butterworth-Bandstop Filter.

Filter data `A` removing data between frequencies `freqmin` to `freqmax` using
`corners` corners.

# Arguments
- `A::AbstractArray`: Data to filter
- `freqmin::Real`: Stop band low corner frequency.
- `freqmax::Real`: Stop band high corner frequency.
- `fs::Real`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function bandstop!(A::AbstractArray,freqmin::Real,freqmax::Real,fs::Real;
                   corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    low = freqmin / fe
    high = freqmax / fe

    # warn if above Nyquist frequency
    if high > 1
        @warn "Selected high corner frequency ($freqmax) is"
        "above Nyquist ($fe). Setting Nyquist as high corner."
        freqmax = fe
    end

    # throw error if low above Nyquist frequency
    if low > 1
        ArgumentError("Selected low corner frequency is above Nyquist.")
    end

    # create filter
    responsetype = Bandstop(freqmin, freqmax; fs=fs)
    designmethod = Butterworth(corners)
    if zerophase
        A[:] = filtfilt(digitalfilter(responsetype, designmethod), A)
    else
        A[:] = filt(digitalfilter(responsetype, designmethod), A)
    end

    return nothing
end
bandstop(A::AbstractArray,freqmin::Real, freqmax::Real; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(A);bandstop!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

function bandstop!(C::SeisChannel, freqmin::Real, freqmax::Real;
                   corners::Int=4, zerophase::Bool=false)
    bandstop!(C.x,freqmin,freqmax,C.fs,corners=corners,zerophase=zerophase)
    return nothing
end
bandstop(C::SeisChannel,freqmin::Real, freqmax::Real; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(C);bandstop!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

function bandstop!(S::SeisData, freqmin::Real, freqmax::Real;
                   corners::Int=4, zerophase::Bool=false)
    @inbounds for i = 1:S.n
        bandstop!(S[i],freqmin,freqmax,corners=corners,zerophase=zerophase)
    end
    return nothing
end

bandstop(S::SeisData,freqmin::Real, freqmax::Real; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(S);bandstop!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

"""
    lowpass(A,freq,fs,corners=4,zerophase=false)

Butterworth-Lowpass Filter.

Filter data `A` over certain frequency `freq` using `corners` corners.

# Arguments
- `A::AbstractArray`: Data to filter
- `freq::Real`: Filter corner frequency.
- `fs::Real`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function lowpass!(A::AbstractArray,freq::Real,fs::Real; corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    f = freq / fe

    # warn if above Nyquist frequency
    if f >= 1
        @warn """Selected corner frequency ($freq) is
        above Nyquist ($fe). Setting Nyquist as high corner."""
        freq = fe - 1. / fs
    end

    # create filter
    responsetype = Lowpass(freq; fs=fs)
    designmethod = Butterworth(corners)
    if zerophase
        A[:] = filtfilt(digitalfilter(responsetype, designmethod), A)
    else
        A[:] = filt(digitalfilter(responsetype, designmethod), A)
    end
    return nothing
end
lowpass(A::AbstractArray,freq::Real; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(A);lowpass!(U,freq,corners=corners,zerophase=zerophase);
        return U)

function lowpass!(C::SeisChannel, freq::Real;corners::Int=4, zerophase::Bool=false)
    lowpass!(C.x,freq,C.fs,corners=corners,zerophase=zerophase)
    return nothing
end
lowpass(C::SeisChannel,freq::Real; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(C);lowpass!(U,freq,corners=corners,zerophase=zerophase);
        return U)

function lowpass!(S::SeisData, freq::Real;corners::Int=4, zerophase::Bool=false)
    @inbounds for i = 1:S.n
        lowpass!(S[i],freq,corners=corners,zerophase=zerophase)
    end
    return nothing
end

lowpass(S::SeisData,freq::Real; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(S);lowpass!(U,freq,corners=corners,zerophase=zerophase);
       return U)

"""
    highpass(A,freq,fs,corners=4,zerophase=false)

Butterworth-Highpass Filter.

Filter data `A` removing data below certain frequency `freq` using `corners` corners.

# Arguments
- `A::AbstractArray`: Data to filter
- `freq::Real`: Filter corner frequency.
- `fs::Real`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function highpass!(A::AbstractArray,freq::Real,fs::Real; corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    f = freq / fe

    # warn if above Nyquist frequency
    if f > 1
        ArgumentError("Selected low corner frequency is above Nyquist.")
    end

    # create filter
    responsetype = Highpass(freq; fs=fs)
    designmethod = Butterworth(corners)
    if zerophase
        A[:] = filtfilt(digitalfilter(responsetype, designmethod), A)
    else
        A[:] = filt(digitalfilter(responsetype, designmethod), A)
    end
    return nothing
end
highpass(A::AbstractArray,freq::Real; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(A);highpass!(U,freq,corners=corners,zerophase=zerophase);
        return U)

function highpass!(C::SeisChannel, freq::Real;corners::Int=4, zerophase::Bool=false)
    highpass!(C.x,freq,C.fs,corners=corners,zerophase=zerophase)
    return nothing
end
highpass(C::SeisChannel,freq::Real; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(C);highpass!(U,freq,corners=corners,zerophase=zerophase);
        return U)

function highpass!(S::SeisData, freq::Real;corners::Int=4, zerophase::Bool=false)
    @inbounds for i = 1:S.n
        highpass!(S[i],freq,corners=corners,zerophase=zerophase)
    end
    return nothing
end
highpass(S::SeisData,freq::Real; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(S);highpass!(U,freq,corners=corners,zerophase=zerophase);
       return U)

"""
   taper!(A,fs; max_percentage=0.05, type="hann", max_length=20)

Taper a time series `A` with sampling_rate `fs`.
Defaults to 'hann' window. Uses smallest of `max_percentage` * `fs`
or `max_length`.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
- `max_percentage::float`: Decimal percentage of taper at one end (ranging
   from 0. to 0.5).
- `max_length::Real`: Length of taper at one end in seconds.
"""
function taper!(A::AbstractArray, fs::Real; max_percentage::Real=0.05,
   max_length::Real=20)
   N = length(A)
   wlen = min(Int(N * max_percentage), Int(max_length * fs), Int(N/2))
   taper_sides = [-hanning(2 * wlen -1, zerophase=true) .+ 1;0]
   A[1:wlen] .= A[1:wlen] .* taper_sides[1:wlen]
   A[end-wlen+1:end] .= A[end-wlen+1:end] .* taper_sides[wlen+1:end]
   return nothing
end
taper(A::AbstractArray, fs::Real; max_percentage::Real=0.05,
   max_length::Real=20) = (U = deepcopy(A);taper!(U,fs,
   max_percentage=max_percentage,max_length=max_length);return U)

function taper!(C::SeisChannel; max_percentage::Real=0.05, max_length::Real=20)
    taper!(C.x,C.fs,max_percentage=max_percentage,max_length=max_length)
    return nothing
end
taper(C::SeisChannel; max_percentage::Real=0.05, max_length::Real=20) =
     (U = deepcopy(A);taper!(U,max_percentage=max_percentage,
      max_length=max_length);return U)

function taper!(S::SeisData; max_percentage::Real=0.05, max_length::Real=20)
    @inbounds for i = 1:S.n
        taper!(S[i],max_percentage=max_percentage, max_length=max_length)
    end
  return nothing
end
taper(S::SeisData; max_percentage::Real=0.05, max_length::Real=20) =
   (U = deepcopy(S);taper!(U,max_percentage=max_percentage,
    max_length=max_length);return U)

"""
    envelope(A)

Envelope of a function.

Computes the upper and lower envelopes of the given function.

# Arguments
- `A::AbstractArray`: Data to make envelope of.
"""
function envelope(A::AbstractArray)
    Amean = mean(A)
    Acentered = A .- Amean
    env = abs.(hilbert(Acentered))
    upper = env .+ Amean
    lower = -env .+ Amean
    return upper, lower
end
