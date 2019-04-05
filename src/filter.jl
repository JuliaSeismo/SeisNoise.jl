module Filter
# filter functions
# Julia translation of obspy.signal.filter
using ..SeisJul
export bandpass, bandstop, lowpass, highpass, envelope

"""
    bandpass(A,freqmin,freqmax,fs,corners=4,zerophase=false)

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
function bandpass(A::AbstractArray, freqmin::Real, freqmax::Real, fs::Real; corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    low = freqmin / fe
    high = freqmax / fe

    # warn if above Nyquist frequency
    if high - oneunit(high) > -1e-6
        @warn "Selected high corner frequency ($freqmax) of bandpass is at or
        above Nyquist ($fe). Applying a high-pass instead."
        return highpass(A,freqmin,fs,corners=corners,zerophase=zerophase)
    end

    # throw error if low above Nyquist frequency
    if low > 1
        ArgumentError("Selected low corner frequency is above Nyquist.")
    end

    # create filter
    responsetype = Bandpass(freqmin, freqmax; fs=fs)
    designmethod = Butterworth(corners)
    if zerophase
        X = filtfilt(digitalfilter(responsetype, designmethod), A)
    else
        X = filt(digitalfilter(responsetype, designmethod), A)
    end

    return X
end

"""
    bandstop(A,freqmin,freqmax,fs,corners=4,zerophase=false)

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
function bandstop(A::AbstractArray,freqmin::Real,freqmax::Real,fs::Real; corners::Int=4, zerophase::Bool=false)
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
        X = filtfilt(digitalfilter(responsetype, designmethod), A)
    else
        X = filt(digitalfilter(responsetype, designmethod), A)
    end

    return X
end

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
function lowpass(A::AbstractArray,freq::Real,fs::Real; corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    f = freq / fe

    # warn if above Nyquist frequency
    if f > 1
        @warn "Selected corner frequency ($freq) is"
        "above Nyquist ($fe). Setting Nyquist as high corner."
        freq = fe
    end

    # create filter
    responsetype = Lowpass(freq; fs=fs)
    designmethod = Butterworth(corners)
    if zerophase
        X = filtfilt(digitalfilter(responsetype, designmethod), A)
    else
        X = filt(digitalfilter(responsetype, designmethod), A)
    end

    return X
end

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
function highpass(A::AbstractArray,freq::Real,fs::Real; corners::Int=4, zerophase::Bool=false)
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
        X = filtfilt(digitalfilter(responsetype, designmethod), A)
    else
        X = filt(digitalfilter(responsetype, designmethod), A)
    end

    return X
end

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

end
