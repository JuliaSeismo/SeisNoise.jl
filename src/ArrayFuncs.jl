export detrend, detrend!, demean, demean!, bandpass, bandpass!, bandstop, bandstop!,
       lowpass, lowpass!, highpass, highpass!, phase!, phase
# Signal processing functions for arrays (rather than SeisData or SeisChannel)

"""
    detrend!(X::AbstractArray{<:Union{Float32,Float64},1})

Remove linear trend from array `X` using least-squares regression.
"""
function detrend!(X::AbstractArray{<:Union{Float32,Float64},1})
    N = length(X)
    A = ones(N,2)
    A[:,1] = Array(1:N) ./ N
    coeff = A \ X
    X[:] .= X .- A *coeff
    return nothing
end
detrend(A::AbstractArray{<:Union{Float32,Float64},1}) = (U = deepcopy(A);
        detrend!(U);return U)

"""
    detrend!(X::AbstractArray{<:Union{Float32,Float64},2})

Remove linear trend from columns of `X` using least-squares regression.
"""
function detrend!(X::AbstractArray{<:Union{Float32,Float64},2})
    M,N = size(X)
    A = ones(M,2)
    A[:,1] = Array(1:M) ./ M

    # solve least-squares through qr decomposition
    Q,R = qr(A)
    rq = inv(factorize(R)) * Q'
    for ii = 1:N
        coeff = rq * X[:,ii]
        X[:,ii] .-=  A *coeff
    end
    return nothing
end
detrend(A::AbstractArray{<:Union{Float32,Float64},2}) = (U = deepcopy(A);
        detrend!(U);return U)

"""
    demean!(A::AbstractArray{<:Union{Float32,Float64},1})

Remove mean from array `A`.
"""
function demean!(A::AbstractArray{<:Union{Float32,Float64},1})
      μ = mean(A)
      for ii = 1:length(A)
        A[ii] -= μ
      end
  return nothing
end
demean(A::AbstractArray{<:Union{Float32,Float64},1}) = (U = deepcopy(A);
       demean!(U);return U)

"""
   demean!(A::AbstractArray{<:Union{Float32,Float64},2})

Remove mean from columns of array `A`.
"""
function demean!(A::AbstractArray{<:Union{Float32,Float64},2})
      M,N = size(A)
      for ii = 1:N
        μ = mean(A[:,ii])
        for jj = 1:M
          A[jj,ii] -= μ
        end
      end
  return nothing
end
demean(A::AbstractArray{<:Union{Float32,Float64},2}) = (U = deepcopy(A);
       demean!(U);return U)

"""
   taper!(A,fs; max_percentage=0.05, max_length=20.)

Taper a time series `A` with sampling_rate `fs`.
Defaults to 'hann' window. Uses smallest of `max_percentage` * `fs`
or `max_length`.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Float64`: Sampling rate of time series `A`.
- `max_percentage::float`: Decimal percentage of taper at one end (ranging
   from 0. to 0.5).
- `max_length::Float64`: Length of taper at one end in seconds.
"""
function taper!(A::AbstractArray{<:Union{Float32,Float64},1}, fs::Float64;
                max_percentage::Float64=0.05, max_length::Float64=20.)
   N = length(A)
   wlen = min(Int(floor(N * max_percentage)), Int(floor(max_length * fs)), Int(
         floor(N/2)))
   taper_sides = [-hanning(2 * wlen -1, zerophase=true) .+ 1;0]
   A[1:wlen] .= A[1:wlen] .* taper_sides[1:wlen]
   A[end-wlen+1:end] .= A[end-wlen+1:end] .* taper_sides[wlen+1:end]
   return nothing
end
taper(A::AbstractArray{<:Union{Float32,Float64},1}, fs::Float64;
      max_percentage::Float64=0.05, max_length::Float64=20.) = (U = deepcopy(A);
      taper!(U,fs,max_percentage=max_percentage,max_length=max_length);return U)

function taper!(A::AbstractArray{<:Union{Float32,Float64},2}, fs::Float64;
                max_percentage::Float64=0.05, max_length::Float64=20.)
   M,N = size(A)
   wlen = min(Int(floor(M * max_percentage)), Int(floor(max_length * fs)), Int(
         floor(M/2)))
   taper_sides = [-hanning(2 * wlen -1, zerophase=true) .+ 1;0]
   for ii = 1:N
       A[1:wlen,ii] .= A[1:wlen,ii] .* taper_sides[1:wlen]
       A[end-wlen+1:end,ii] .= A[end-wlen+1:end,ii] .* taper_sides[wlen+1:end]
   end
   return nothing
end
taper(A::AbstractArray{<:Union{Float32,Float64},2}, fs::Float64;
       max_percentage::Float64=0.05,max_length::Float64=20.) = (U = deepcopy(A);
       taper!(U,fs,max_percentage=max_percentage,max_length=max_length);return U)

"""
    phase!(A::AbstractArray)

Extract instantaneous phase from signal A.

For time series `A`, its analytic representation ``S = A + H(A)``, where
``H(A)`` is the Hilbert transform of `A`. The instantaneous phase ``e^{iθ}``
of `A` is given by dividing ``S`` by its modulus: ``e^{iθ} = \\frac{S}{|S|}``
For more information on Phase Cross-Correlation, see:
[Ventosa et al., 2019](https://pubs.geoscienceworld.org/ssa/srl/article-standard/570273/towards-the-processing-of-large-data-volumes-with).
"""
function phase!(A::AbstractArray)
    A .= angle.(hilbert(A))
    return nothing
end
phase(A::AbstractArray) = (U = deepcopy(A);phase!(U);return U)


"""
   bandpass!(A,freqmin,freqmax,fs,corners=4,zerophase=false)

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
function bandpass!(A::AbstractArray{<:Union{Float32,Float64},1},
                   freqmin::Float64, freqmax::Float64, fs::Float64;
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
bandpass(A::AbstractArray{<:Union{Float32,Float64},1},freqmin::Float64,
         freqmax::Float64, fs::Float64; corners::Int=4,zerophase::Bool=false) =
         (U = deepcopy(A);bandpass!(U,freqmin,freqmax, fs, corners=corners,
         zerophase=zerophase);return U)

function bandpass!(A::AbstractArray{<:Union{Float32,Float64},2},
                   freqmin::Float64, freqmax::Float64, fs::Float64;
                   corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    low = freqmin / fe
    high = freqmax / fe
    M,N = size(A)

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
        for ii = 1:N
            A[:,ii] = filtfilt(digitalfilter(responsetype, designmethod), A[:,ii])
        end
    else
        for ii = 1:N
            A[:,ii] = filt(digitalfilter(responsetype, designmethod), A[:,ii])
        end
    end

    return nothing
end
bandpass(A::AbstractArray{<:Union{Float32,Float64},2},freqmin::Float64,
         freqmax::Float64, fs::Float64; corners::Int=4,zerophase::Bool=false) =
         (U = deepcopy(A);bandpass!(U,freqmin,freqmax,fs,corners=corners,
         zerophase=zerophase);return U)


 """
     bandstop!(A,freqmin,freqmax,fs,corners=4,zerophase=false)

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
function bandstop!(A::AbstractArray{<:Union{Float32,Float64},1},
                    freqmin::Float64,freqmax::Float64,fs::Float64;
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
bandstop(A::AbstractArray{<:Union{Float32,Float64},1},freqmin::Float64,
         freqmax::Float64, fs::Float64; corners::Int=4,zerophase::Bool=false) =
         (U = deepcopy(A);bandstop!(U,freqmin,freqmax,corners=corners,
         zerophase=zerophase);return U)

function bandstop!(A::AbstractArray{<:Union{Float32,Float64},2},
                   freqmin::Float64,freqmax::Float64,fs::Float64;
                   corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    low = freqmin / fe
    high = freqmax / fe
    M,N = size(A)

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
        for ii = 1:N
            A[:,ii] = filtfilt(digitalfilter(responsetype, designmethod), A[:,ii])
        end
    else
        for ii = 1:N
            A[:,ii] = filt(digitalfilter(responsetype, designmethod), A[:,ii])
        end
    end

    return nothing
end
bandstop(A::AbstractArray{<:Union{Float32,Float64},2},freqmin::Float64,
         freqmax::Float64, fs::Float64; corners::Int=4,zerophase::Bool=false) =
         (U = deepcopy(A);bandstop!(U,freqmin,freqmax,corners=corners,
         zerophase=zerophase);return U)

"""
lowpass(A,freq,fs,corners=4,zerophase=false)

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
function lowpass!(A::AbstractArray{<:Union{Float32,Float64},1},freq::Float64,
                  fs::Float64; corners::Int=4, zerophase::Bool=false)
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
lowpass(A::AbstractArray{<:Union{Float32,Float64},1},freq::Float64, fs::Float64;
        corners::Int=4,zerophase::Bool=false) = (U = deepcopy(A);
        lowpass!(U,freq,fs,corners=corners,zerophase=zerophase);return U)

function lowpass!(A::AbstractArray{<:Union{Float32,Float64},2},freq::Float64,
                  fs::Float64; corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    f = freq / fe
    M,N = size(A)

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
        for ii = 1:N
            A[:,ii] = filtfilt(digitalfilter(responsetype, designmethod), A[:,ii])
        end
    else
        for ii = 1:N
            A[:,ii] = filt(digitalfilter(responsetype, designmethod), A[:,ii])
        end
    end
    return nothing
end
lowpass(A::AbstractArray{<:Union{Float32,Float64},2},freq::Float64, fs::Float64;
        corners::Int=4,zerophase::Bool=false) = (U = deepcopy(A);
        lowpass!(U,freq,fs,corners=corners,zerophase=zerophase);return U)

"""
    highpass(A,freq,fs,corners=4,zerophase=false)

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
function highpass!(A::AbstractArray{<:Union{Float32,Float64},1},freq::Float64,
                   fs::Float64; corners::Int=4, zerophase::Bool=false)
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
highpass(A::AbstractArray{<:Union{Float32,Float64},1},freq::Float64,fs::Float64;
         corners::Int=4,zerophase::Bool=false) = (U = deepcopy(A);
         highpass!(U,freq,fs,corners=corners,zerophase=zerophase);return U)

function highpass!(A::AbstractArray{<:Union{Float32,Float64},2},freq::Float64,
        fs::Float64; corners::Int=4, zerophase::Bool=false)
    fe = 0.5 * fs
    f = freq / fe
    M,N = size(A)

    # warn if above Nyquist frequency
    if f > 1
        ArgumentError("Selected low corner frequency is above Nyquist.")
    end

    # create filter
    responsetype = Highpass(freq; fs=fs)
    designmethod = Butterworth(corners)
    if zerophase
        for ii = 1:N
            A[:,ii] = filtfilt(digitalfilter(responsetype, designmethod), A[:,ii])
        end
    else
        for ii = 1:N
            A[:,ii] = filt(digitalfilter(responsetype, designmethod), A[:,ii])
        end
    end
    return nothing
end
highpass(A::AbstractArray{<:Union{Float32,Float64},2},freq::Float64,fs::Float64;
corners::Int=4,zerophase::Bool=false) = (U = deepcopy(A);
highpass!(U,freq,fs,corners=corners,zerophase=zerophase);return U)
