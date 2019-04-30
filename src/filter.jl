export bandpass, bandpass!, bandstop, bandstop!, lowpass, lowpass!,
        highpass, highpass!, taper, taper!, envelope

"""
    bandpass!(C,freqmin,freqmax,fs,corners=4,zerophase=false)

Butterworth-Bandpass Filter.

Filter data in `C` from `freqmin` to `freqmax` using `corners` corners.

# Arguments
- `C::SeisChannel`: SeisChannel to filter.
- `freqmin::Float64`: Pass band low corner frequency.
- `freqmax::Float64`: Pass band high corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function bandpass!(C::SeisChannel, freqmin::Float64, freqmax::Float64;
                   corners::Int=4, zerophase::Bool=false)
    bandpass!(C.x,freqmin,freqmax,C.fs,corners=corners,zerophase=zerophase)
    return nothing
end
bandpass(C::SeisChannel,freqmin::Float64, freqmax::Float64; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(C);bandpass!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

"""
  bandpass!(S,freqmin,freqmax,fs,corners=4,zerophase=false)

Butterworth-Bandpass Filter.

Filter channels in `S` from `freqmin` to `freqmax` using `corners` corners.

# Arguments
- `S::SeisData`: SeisData to filter.
- `freqmin::Float64`: Pass band low corner frequency.
- `freqmax::Float64`: Pass band high corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function bandpass!(S::SeisData, freqmin::Float64, freqmax::Float64;
                   corners::Int=4, zerophase::Bool=false)
    @inbounds for i = 1:S.n
        bandpass!(S[i].x,freqmin,freqmax,S[i].fs,corners=corners,
                             zerophase=zerophase)
    end
    return nothing
end
bandpass(S::SeisData,freqmin::Float64, freqmax::Float64; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(S);bandpass!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

"""
    bandstop!(C,freqmin,freqmax,fs,corners=4,zerophase=false)

Butterworth-Bandstop Filter.

Filter data in `C` removing data between frequencies `freqmin` to `freqmax` using
`corners` corners.

# Arguments
- `C::SeisChannel`: SeisChannel to filter.
- `freqmin::Float64`: Stop band low corner frequency.
- `freqmax::Float64`: Stop band high corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function bandstop!(C::SeisChannel, freqmin::Float64, freqmax::Float64;
                   corners::Int=4, zerophase::Bool=false)
    bandstop!(C.x,freqmin,freqmax,C.fs,corners=corners,zerophase=zerophase)
    return nothing
end
bandstop(C::SeisChannel,freqmin::Float64, freqmax::Float64; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(C);bandstop!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

function bandstop!(S::SeisData, freqmin::Float64, freqmax::Float64;
                   corners::Int=4, zerophase::Bool=false)
    @inbounds for i = 1:S.n
        bandstop!(S[i].x,freqmin,freqmax,s[i].fs,corners=corners,zerophase=zerophase)
    end
    return nothing
end

"""
    bandstop!(S,freqmin,freqmax,fs,corners=4,zerophase=false)

Butterworth-Bandstop Filter.

Filter channels in `S` removing data between frequencies `freqmin` to `freqmax` using
`corners` corners.

# Arguments
- `S::SeisData`: SeisData to filter.
- `freqmin::Float64`: Stop band low corner frequency.
- `freqmax::Float64`: Stop band high corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
bandstop(S::SeisData,freqmin::Float64, freqmax::Float64; corners::Int=4,
         zerophase::Bool=false) = (U = deepcopy(S);bandstop!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)

"""
    lowpass(C,freq,fs,corners=4,zerophase=false)

Butterworth-Lowpass Filter.

Filter data in `C` over certain frequency `freq` using `corners` corners.

# Arguments
- `C::SeisChannel`: SeisChannel to filter.
- `freq::Float64`: Filter corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function lowpass!(C::SeisChannel, freq::Float64;corners::Int=4, zerophase::Bool=false)
    lowpass!(C.x,freq,C.fs,corners=corners,zerophase=zerophase)
    return nothing
end
lowpass(C::SeisChannel,freq::Float64; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(C);lowpass!(U,freq,corners=corners,zerophase=zerophase);
        return U)

"""
    lowpass(S,freq,fs,corners=4,zerophase=false)

Butterworth-Lowpass Filter.

Filter channels in `S` over certain frequency `freq` using `corners` corners.

# Arguments
- `S::SeisData`: SeisData to filter.
- `freq::Float64`: Filter corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function lowpass!(S::SeisData, freq::Float64;corners::Int=4, zerophase::Bool=false)
    @inbounds for i = 1:S.n
        lowpass!(S[i].x,freq,S[i].fs,corners=corners,
                            zerophase=zerophase)
    end
    return nothing
end

lowpass(S::SeisData,freq::Float64; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(S);lowpass!(U,freq,corners=corners,zerophase=zerophase);
       return U)

"""
    highpass(C,freq,fs,corners=4,zerophase=false)

Butterworth-Highpass Filter.

Filter data in `C` removing data below certain frequency `freq` using `corners` corners.

# Arguments
- `C::SeisChannel`: SeisChannel to filter.
- `freq::Float64`: Filter corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function highpass!(C::SeisChannel, freq::Float64;corners::Int=4, zerophase::Bool=false)
    highpass!(C.x,freq,C.fs,corners=corners,zerophase=zerophase)
    return nothing
end
highpass(C::SeisChannel,freq::Float64; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(C);highpass!(U,freq,corners=corners,zerophase=zerophase);
        return U)

"""
    highpass(S,freq,fs,corners=4,zerophase=false)

Butterworth-Highpass Filter.

Filter channels in `S` removing data below certain frequency `freq` using `corners` corners.

# Arguments
- `S::SeisData`: SeisData to filter.
- `freq::Float64`: Filter corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function highpass!(S::SeisData, freq::Float64;corners::Int=4, zerophase::Bool=false)
    @inbounds for i = 1:S.n
        highpass!(S[i].x,freq,S[i].fs,corners=corners,
                             zerophase=zerophase)
    end
    return nothing
end
highpass(S::SeisData,freq::Float64; corners::Int=4,zerophase::Bool=false) =
       (U = deepcopy(S);highpass!(U,freq,corners=corners,zerophase=zerophase);
       return U)

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
