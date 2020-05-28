export bandpass, bandpass!, bandstop, bandstop!, lowpass, lowpass!,
        highpass, highpass!, taper, taper!, envelope

"""
   bandpass!(A,freqmin,freqmax,fs,corners=4,zerophase=true)

Butterworth-Bandpass Filter.

Filter data `A` from `freqmin` to `freqmax` using `corners` corners.

# Arguments
- `A::AbstractArray`: Data to filter
- `freqmin::Float64`: Pass band low corner frequency.
- `freqmax::Float64`: Pass band high corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function bandpass!(A::AbstractArray{<:AbstractFloat},
                   freqmin::Float64, freqmax::Float64, fs::Float64;
                   corners::Int=4, zerophase::Bool=true)
   T = eltype(A)
   fe = T(0.5 * fs)
   low = T(freqmin / fe)
   high = T(freqmax / fe)

   # warn if above Nyquist frequency
   if high - oneunit(high) > -1e-6
       @warn "Selected high corner frequency ($freqmax) of bandpass is at or
       above Nyquist ($fe). Applying a high-pass instead."
       highpass!(A,freqmin,fs,corners=corners,zerophase=zerophase)
       return nothing
   end

   # throw error if low above Nyquist frequency
   if low >= 1
       throw(ArgumentError("Selected low corner frequency is above Nyquist."))
   end

   # create filter
   responsetype = Bandpass(T(freqmin), T(freqmax); fs=fs)
   designmethod = Butterworth(T,corners)

   # use gpu-specific kernel if on the GPU
   if isa(A,AbstractGPUArray)
       gpufilter!(A,responsetype,designmethod)
       return nothing
   end

   # filter if on the CPU
   if zerophase
       A[:,:] .= filtfilt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
   else
       A[:,:] .= filt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
   end

   return nothing
end
bandpass(A::AbstractArray{<:AbstractFloat},freqmin::Float64,
         freqmax::Float64, fs::Float64; corners::Int=4,zerophase::Bool=true) =
         (U = deepcopy(A);bandpass!(U,freqmin,freqmax, fs, corners=corners,
          zerophase=zerophase);return U)
bandpass!(R::RawData,freqmin::Float64,freqmax::Float64;
          corners::Int=4,zerophase::Bool=true) = (bandpass!(R.x,freqmin,freqmax,
          R.fs,corners=corners,zerophase=zerophase);setfield!(R,:freqmin,freqmin);
          setfield!(R,:freqmax,min(freqmax,R.fs/2));return nothing)
bandpass(R::RawData,freqmin::Float64,freqmax::Float64;
        corners::Int=4,zerophase::Bool=true) = (U = deepcopy(R);bandpass!(U.x,
        freqmin,freqmax,U.fs,corners=corners,zerophase=zerophase);
        setfield!(U,:freqmin,freqmin);
        setfield!(U,:freqmax,min(freqmax,U.fs/2));return U)
bandpass!(C::CorrData,freqmin::Float64,freqmax::Float64;
          corners::Int=4,zerophase::Bool=true) = (bandpass!(C.corr,freqmin,freqmax,
          C.fs,corners=corners,zerophase=zerophase);setfield!(C,:freqmin,freqmin);
          setfield!(C,:freqmax,min(freqmax,C.fs/2));return nothing)
bandpass(C::CorrData,freqmin::Float64,freqmax::Float64;
        corners::Int=4,zerophase::Bool=true) = (U = deepcopy(C);bandpass!(U.corr,
        freqmin,freqmax,U.fs,corners=corners,zerophase=zerophase);
        setfield!(U,:freqmin,freqmin);
        setfield!(U,:freqmax,min(freqmax,U.fs/2));return U)
bandpass!(C::SeisChannel, freqmin::Float64, freqmax::Float64;
                   corners::Int=4, zerophase::Bool=true) = filtfilt!(C,
                   fl=freqmin,fh=freqmax,np=corners,rt="Bandpass")
bandpass(C::SeisChannel,freqmin::Float64, freqmax::Float64; corners::Int=4,
         zerophase::Bool=true) = (U = deepcopy(C);bandpass!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)
bandpass!(S::SeisData, freqmin::Float64, freqmax::Float64;
                   corners::Int=4, zerophase::Bool=true) = filtfilt!(S,
                   fl=freqmin,fh=freqmax,np=corners,rt="Bandpass")
bandpass(S::SeisData, freqmin::Float64, freqmax::Float64; corners::Int=4,
         zerophase::Bool=true) = filtfilt(S,fl=freqmin,fh=freqmax,np=corners,rt="Bandpass")

"""
  bandstop!(A,freqmin,freqmax,fs,corners=4,zerophase=true)

Butterworth-Bandstop Filter.

Filter data `A` removing data between frequencies `freqmin` to `freqmax` using
`corners` corners.

# Arguments
- `A::AbstractArray`: Data to filter
- `freqmin::Float64`: Stop band low corner frequency.
- `freqmax::Float64`: Stop band high corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function bandstop!(A::AbstractArray{<:AbstractFloat},
                 freqmin::Float64,freqmax::Float64,fs::Float64;
                 corners::Int=4, zerophase::Bool=true)
    T = eltype(A)
    fe = T(0.5 * fs)
    low = T(freqmin / fe)
    high = T(freqmax / fe)

    # throw error if high above Nyquist frequency
    if high >= 1
        throw(ArgumentError("Selected high corner frequency $freqmax Hz is above Nyquist $(fs/2) Hz."))
    end

    # throw error if low above Nyquist frequency
    if low >= 1
        throw(ArgumentError("Selected low corner frequency $freqmin Hz is above Nyquist $(fs/2) Hz."))
    end

    # create filter
    responsetype = Bandstop(T(freqmin), T(freqmax); fs=fs)
    designmethod = Butterworth(T,corners)

    # use gpu-specific kernel if on the GPU
    if isa(A,AbstractGPUArray)
        gpufilter!(A,responsetype,designmethod)
        return nothing
    end

    # filter if on the CPU
    if zerophase
        A[:,:] .= filtfilt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
    else
        A[:,:] .= filt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
    end

    return nothing
end
bandstop(A::AbstractArray{<:AbstractFloat},freqmin::Float64,
      freqmax::Float64, fs::Float64; corners::Int=4,zerophase::Bool=true) =
      (U = deepcopy(A);bandstop!(U,freqmin,freqmax,fs,corners=corners,
      zerophase=zerophase);return U)
bandstop!(R::RawData,freqmin::Float64,freqmax::Float64;
 corners::Int=4,zerophase::Bool=true) = bandstop!(R.x,freqmin,freqmax,
 R.fs,corners=corners,zerophase=zerophase)
bandstop(R::RawData,freqmin::Float64,freqmax::Float64;
corners::Int=4,zerophase::Bool=true) = (U = deepcopy(R);bandstop!(U.x,
freqmin,freqmax,U.fs,corners=corners,zerophase=zerophase);return U)
bandstop!(C::CorrData,freqmin::Float64,freqmax::Float64;
 corners::Int=4,zerophase::Bool=true) = bandstop!(C.corr,freqmin,freqmax,
 C.fs,corners=corners,zerophase=zerophase)
bandstop(C::CorrData,freqmin::Float64,freqmax::Float64;
corners::Int=4,zerophase::Bool=true) = (U = deepcopy(C);bandstop!(U.corr,
freqmin,freqmax,U.fs,corners=corners,zerophase=zerophase);return U)
bandstop!(C::SeisChannel, freqmin::Float64, freqmax::Float64;
                   corners::Int=4, zerophase::Bool=true) = filtfilt!(C,
                   fl=freqmin,fh=freqmax,np=corners,rt="Bandstop")
bandstop(C::SeisChannel,freqmin::Float64, freqmax::Float64; corners::Int=4,
         zerophase::Bool=true) = filtfilt(C,
         fl=freqmin,fh=freqmax,np=corners,rt="Bandstop")
bandstop!(S::SeisData, freqmin::Float64, freqmax::Float64;
                   corners::Int=4, zerophase::Bool=true) = filtfilt!(S,
                   fl=freqmin,fh=freqmax,np=corners,rt="Bandstop")
bandstop(S::SeisData,freqmin::Float64, freqmax::Float64; corners::Int=4,
         zerophase::Bool=true) = filtfilt(S,
         fl=freqmin,fh=freqmax,np=corners,rt="Bandstop")

"""
lowpass(A,freq,fs,corners=4,zerophase=true)

Butterworth-Lowpass Filter.

Filter data `A` over certain frequency `freq` using `corners` corners.

# Arguments
- `A::AbstractArray`: Data to filter
- `freq::Float64`: Filter corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function lowpass!(A::AbstractArray{<:AbstractFloat},freq::Float64,
                fs::Float64; corners::Int=4, zerophase::Bool=true)
    fe = 0.5 * fs
    f = freq / fe
    T = eltype(A)

    # warn if above Nyquist frequency
    if f >= 1
        @warn """Selected corner frequency ($freq) is
        above Nyquist ($fe). No filter applied."""
        return nothing
    end

    # create filter
    responsetype = Lowpass(T(freq); fs=fs)
    designmethod = Butterworth(T,corners)

    # use gpu-specific kernel if on the GPU
    if isa(A,AbstractGPUArray)
        gpufilter!(A,responsetype,designmethod)
        return nothing
    end

    # filter if on the CPU
    if zerophase
        A[:,:] .= filtfilt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
    else
        A[:,:] .= filt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
    end
    return nothing
end
lowpass(A::AbstractArray{<:AbstractFloat},freq::Float64, fs::Float64;
      corners::Int=4,zerophase::Bool=true) = (U = deepcopy(A);
      lowpass!(U,freq,fs,corners=corners,zerophase=zerophase);return U)
lowpass!(R::RawData,freq::Float64; corners::Int=4,
       zerophase::Bool=true) = (lowpass!(R.x,freq,R.fs,corners=corners,
       zerophase=zerophase);setfield!(R,:freqmax,min(freq,R.fs/2));
       return nothing)
lowpass(R::RawData,freq::Float64; corners::Int=4,
       zerophase::Bool=true) = (U = deepcopy(R);lowpass!(U.x,freq,U.fs,
       corners=corners,zerophase=zerophase);
       setfield!(U,:freqmax,min(freq,U.fs/2));return U)
lowpass!(C::CorrData,freq::Float64; corners::Int=4,
      zerophase::Bool=true) = (lowpass!(C.corr,freq,C.fs,corners=corners,
      zerophase=zerophase);setfield!(C,:freqmax,min(freq,C.fs/2));
      return nothing)
lowpass(C::CorrData,freq::Float64; corners::Int=4,
      zerophase::Bool=true) = (U = deepcopy(C);lowpass!(U.corr,freq,U.fs,
      corners=corners,zerophase=zerophase);
      setfield!(U,:freqmax,min(freq,U.fs/2));return U)
lowpass!(C::SeisChannel, freq::Float64;corners::Int=4, zerophase::Bool=true) =
    filtfilt!(C,fh=freq,np=corners,rt="Lowpass")
lowpass(C::SeisChannel,freq::Float64; corners::Int=4,zerophase::Bool=true) =
    filtfilt(C,fh=freq,np=corners,rt="Lowpass")
lowpass!(S::SeisData, freq::Float64;corners::Int=4, zerophase::Bool=true) =
    filtfilt!(S,fh=freq,np=corners,rt="Lowpass")
lowpass(S::SeisData,freq::Float64; corners::Int=4,zerophase::Bool=true) =
    filtfilt(S,fh=freq,np=corners,rt="Lowpass")

"""
   highpass(A,freq,fs,corners=4,zerophase=true)

Butterworth-Highpass Filter.

Filter data `A` removing data below certain frequency `freq` using `corners` corners.

# Arguments
- `A::AbstractArray`: Data to filter
- `freq::Float64`: Filter corner frequency.
- `fs::Float64`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.
"""
function highpass!(A::AbstractArray{<:AbstractFloat},freq::Float64,
                  fs::Float64; corners::Int=4, zerophase::Bool=true)
    fe = 0.5 * fs
    f = freq / fe
    T = eltype(A)

    # warn if above Nyquist frequency
    if f >= 1
        throw(ArgumentError("Selected low corner frequency $freq Hz is above Nyquist $(fs/2) Hz."))
    end

    # create filter
    responsetype = Highpass(T(freq); fs=fs)
    designmethod = Butterworth(T,corners)

    # use gpu-specific kernel if on the GPU
    if isa(A,AbstractGPUArray)
        gpufilter!(A,responsetype,designmethod)
        return nothing
    end

    # filter if on the CPU
    if zerophase
        A[:,:] .= filtfilt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
    else
        A[:,:] .= filt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
    end
    return nothing
end
highpass(A::AbstractArray{<:AbstractFloat},freq::Float64,fs::Float64;
        corners::Int=4,zerophase::Bool=true) = (U = deepcopy(A);
        highpass!(U,freq,fs,corners=corners,zerophase=zerophase);return U)
highpass!(R::RawData,freq::Float64; corners::Int=4,
        zerophase::Bool=true) = (highpass!(R.x,freq,R.fs,corners=corners,
        zerophase=zerophase);setfield!(R,:freqmin,freq);return nothing)
highpass(R::RawData,freq::Float64; corners::Int=4,
        zerophase::Bool=true) = (U = deepcopy(R);highpass!(U.x,freq,U.fs,
        corners=corners,zerophase=zerophase);setfield!(U,:freqmin,freq);
        return U)
highpass!(C::CorrData,freq::Float64; corners::Int=4,
        zerophase::Bool=true) = (highpass!(C.corr,freq,C.fs,corners=corners,
        zerophase=zerophase);setfield!(C,:freqmin,freq);return nothing)
highpass(C::CorrData,freq::Float64; corners::Int=4,
        zerophase::Bool=true) = (U = deepcopy(C);highpass!(U.corr,freq,U.fs,
        corners=corners,zerophase=zerophase);setfield!(U,:freqmin,freq);
        return U)
highpass!(C::SeisChannel, freq::Float64;corners::Int=4, zerophase::Bool=true) =
    filtfilt!(C,fl=freq,np=corners,rt="Highpass")
highpass(C::SeisChannel,freq::Float64; corners::Int=4,zerophase::Bool=true) =
    filtfilt(C,fl=freq,np=corners,rt="Highpass")
highpass!(S::SeisData, freq::Float64;corners::Int=4, zerophase::Bool=true) =
    filtfilt!(S,fl=freq,np=corners,rt="Highpass")
highpass(S::SeisData,freq::Float64; corners::Int=4,zerophase::Bool=true) =
    filtfilt(S,fl=freq,np=corners,rt="Highpass")

"""
    envelope(A)

Envelope of a function.

Computes the upper and lower envelopes of the given function.

# Arguments
- `A::AbstractArray`: Data to make envelope of.
"""
function envelope(A::AbstractArray)
    Amean = mean(A,dims=1)
    Acentered = A .- Amean
    env = abs.(hilbert(Acentered))
    upper = env .+ Amean
    lower = -env .+ Amean
    return upper, lower
end


"""
    gpufilter!(A,responsetype,designmethod)

Apply filter to array `A` on the GPU.

# Arguments
- `A::AbstractGPUArray`: Array on the GPU to filter
- `responsetype::FilterType`: DSP.jl filter representation
- `designmethod::ZeroPoleGain`: Filter representation in terms of zeros `z`, poles `p`, and
gain `k`.
"""
function gpufilter!(A::AbstractGPUArray,responsetype::FilterType,designmethod::ZeroPoleGain)
    T = Float64
    N = size(A,1)
    tf = convert(PolynomialRatio, digitalfilter(responsetype, designmethod))

    # get filter coefficients in Float64 to avoid divide by zero
    b = coefb(tf)
    a = coefa(tf)
    newb = CuArrays.zeros(T,N)
    newa = CuArrays.zeros(T,N)
    copyto!(newb,b)
    copyto!(newa,a)
    bafft = rfft(newb) ./ rfft(newa)

    # convert back to Float32 for speed
    bafft = map(x->convert(ComplexF32,x), bafft)

    # apply filter
    A .= irfft(rfft(A,1) .* bafft, N, 1)
    reverse!(A,dims=1)
    A .= irfft(rfft(A,1) .* bafft, N, 1)
    reverse!(A,dims=1)
    return nothing
end
