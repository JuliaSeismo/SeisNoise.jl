# cross-correlation module
export clean_up!, clean_up, correlate, map_cc, corrmap, phasecorrelate
export coherence!, coherence, deconvolution!, deconvolution, whiten, whiten!

"""
    clean_up!(A,freqmin,freqmax,fs)

Demean, detrend, taper and filter time series.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Float64`: Sampling rate of time series `A` in Hz.
- `freqmin::Float64`: Pass band low corner frequency in Hz.
- `freqmax::Float64`: Pass band high corner frequency in Hz.
"""
function clean_up!(A::AbstractArray, freqmin::Real, freqmax::Real, fs::Real;
                   corners::Int=4,zerophase::Bool=true,max_length::Float64=20.)
    detrend!(A)
    taper!(A,fs,max_length=max_length)
    bandpass!(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
    return nothing
end
clean_up(A::AbstractArray, freqmin::Real, freqmax::Real, fs::Real;
         corners::Int=4, zerophase::Bool=true,max_length::Float64=20.) =
         (U = deepcopy(A); clean_up!(U,freqmin,freqmax, fs,
         corners=corners, zerophase=zerophase, max_length=max_length); return U)

clean_up!(C::CorrData,freqmin::Float64,freqmax::Float64; corners::Int=4,
          zerophase::Bool=true,max_length::Float64=20.) = (clean_up!(C.corr,
          freqmin,freqmax,C.fs,corners=corners,zerophase=zerophase,
          max_length=max_length);C.freqmin=freqmin;C.freqmax=freqmax;return nothing)

clean_up(C::CorrData,freqmin::Float64,freqmax::Float64; corners::Int=4,
         zerophase::Bool=true,max_length::Float64=20.) = (U = deepcopy(C);
         clean_up!(U,freqmin,freqmax,corners=corners,zerophase=zerophase,
         max_length=max_length);return U)

clean_up!(R::RawData,freqmin::Float64,freqmax::Float64; corners::Int=4,
          zerophase::Bool=true,max_length::Float64=20.) = (clean_up!(R.x,freqmin,
          freqmax,R.fs,corners=corners,zerophase=true, max_length=max_length);
          R.freqmin=freqmin;R.freqmax=freqmax;return nothing)

clean_up(R::RawData,freqmin::Float64,freqmax::Float64; corners::Int=4,
       zerophase::Bool=true,max_length::Float64=20.) = (U = deepcopy(R);
       clean_up!(U,freqmin,freqmax,corners=corners,zerophase=zerophase,
       max_length=max_length);return U)

"""
    correlate(FFT1, FFT2, N, maxlag, corr_type='cross-correlation')

Cross-correlate ambient noise data in the frequency domain.

Cross-correlation can be done using one of three options:

- Cross-correlation: ``C_{AB}(ω) = u_A(ω) u^∗_B(ω)``
- Coherence: ``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_A(ω) ∣  ∣ u_B(ω) ∣}``
- Deconvolution: ``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_B(ω) ∣^2}``

Smoothing of FFTs for coherence and deconvolution should be done before
cross-correlating.

# Arguments
- `FFT1::AbstractArray`: Complex Array of fourier transform of ambient noise data.
- `FFT2::AbstractArray`: Complex Array of fourier transform of ambient noise data.
- `N::Int`: Number of input data points in time domain, equal to `cc_len` * `fs`.
- `maxlag::Int`: Number of data points in cross-correlation to save,
                 e.g. `maxlag = 2000` will save lag times = -2000/fs:2000/fs s.
- `corr_type::String`: Type of correlation: `cross-correlation`, `coherence` or
                       `deconv`.
"""
function correlate(FFT1::AbstractArray{Complex{T}}, FFT2::AbstractArray{Complex{T}},
                      N::Int, maxlag::Int) where T <: AbstractFloat
    # take inverse fft
    corrT = irfft(conj.(FFT1) .* FFT2,N,1)

    # return corr[-maxlag:maxlag]
    t = vcat(0:Int(N  / 2)-1, -Int(N  / 2):-1)
    ind = findall(abs.(t) .<= maxlag)
    newind = fftshift(ind,1)
    return corrT[newind,:]
end

"""
    phasecorrelate(FFT1, FFT2, maxlag)

Phase Cross-correlate (PCC) ambient noise data in the frequency domain.

# Arguments
- `FFT1::AbstractArray`: Complex Array of fourier transform of ambient noise data.
- `FFT2::AbstractArray`: Complex Array of fourier transform of ambient noise data.
- `N::Int`: Number of input data points in time domain, equal to `cc_len` * `fs`.
- `maxlag::Int`: Number of data points in cross-correlation to save,
                 e.g. `maxlag = 2000` will save lag times = -2000/fs:2000/fs s.
- `corr_type::String`: Type of correlation: `cross-correlation`, `coherence` or
                       `deconv`.
"""
function phasecorrelate(FFT1::AbstractArray{Complex{T}}, FFT2::AbstractArray{Complex{T}},
                      N::Int, maxlag::Int) where T <: AbstractFloat
    # take inverse fft
    corrT = real.(ifft(conj.(FFT1) .* FFT2,1))

    # return corr[-maxlag:maxlag]
    t = vcat(0:Int(N  / 2)-1, -Int(N  / 2):-1)
    ind = findall(abs.(t) .<= maxlag)
    newind = fftshift(ind,1)
    return corrT[newind,:]
end

"""
    correlate(FFT1, FFT2, maxlag,corr_type="CC")

Cross-correlate ambient noise data in the frequency domain.

Cross-correlation can be done using one of two options:

- CC: Cross-correlation, i.e. ``C_{AB}(ω) = u_A(ω) u^∗_B(ω)``
- PCC: Phase cross-correlation, see [Ventosa et al., 2019]

# Arguments
- `FFT1::FFTData`: FFTData object of fft'd ambient noise data.
- `FFT2::FFTData`: FFTData object of fft'd ambient noise data.
- `maxlag::Float64`: Maximum lag time (in seconds) in cross-correlation to save,
                     e.g. `maxlag = 20.` will save lag times = -20.:20. s.
- `corr_type::String`: Type of correlation: `CC` or `PCC`.
"""
function correlate(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;corr_type::String="CC")

    comp = FFT1.name[end] * FFT2.name[end]
    # get intersect of dates; return nothing if no intersect
    inter = intersect(FFT1.t,FFT2.t)
    if length(inter) == 0
        throw(ArgumentError("No common windows for $(FFT1.name)-$(FFT2.name) $(FFT1.id)"))
    end

    ind1 = findall(x -> x ∈ inter, FFT1.t)
    ind2 = findall(x -> x ∈ inter, FFT2.t)
    N = convert(Int,round(FFT1.cc_len * FFT1.fs)) # number of data points
    if uppercase(corr_type) == "CC"
        corr = correlate(@views(FFT1.fft[:,ind1]), @views(FFT2.fft[:,ind2]),
                     N,convert(Int,round(maxlag * FFT1.fs)))
    elseif uppercase(corr_type) == "PCC"
        corr = phasecorrelate(@views(FFT1.fft[:,ind1]), @views(FFT2.fft[:,ind2]),
                     N,convert(Int,round(maxlag * FFT1.fs)))
    else
        throw(ArgumentError("Unrecognized cross-correlation type $corr_type. Options are CC and PCC."))
    end

    rotated = false

    return CorrData(FFT1, FFT2, comp, rotated, corr_type,maxlag, inter, corr)
end

"""
   whiten!(A, freqmin, freqmax, fs, pad=50)

Whiten spectrum of rfft `A` between frequencies `freqmin` and `freqmax`.
Returns the whitened rfft of the time series.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
- `N::Int`: Number of input time domain samples for each rfft.
- `pad::Int`: Number of tapering points outside whitening band.
"""
function whiten!(A::AbstractArray{Complex{Float32}}, freqmin::Float64,
                 freqmax::Float64, fs::Float64,N::Int;pad::Int=50)
   T = real(eltype(A))
   Nrows,Ncols = size(A)

   # get whitening frequencies
   freqvec = FFTW.rfftfreq(N,fs)
   left = findfirst(x -> x >= freqmin, freqvec)
   right = findfirst(freqmax .<= freqvec)
   low, high = left - pad, right + pad

   if low <= 1
       low = 1
       left = low + pad
   end

   if high > length(freqvec)
       high = length(freqvec)- 1
       right = high - pad
   end

   compzero = complex(T(0))
   padarr = similar(A,T,pad)
   padarr .= T(0.):T(pad-1)
   # left zero cut-off
   A[1:low-1,:] .= compzero
   # left tapering
   A[low:left-1,:] .= cos.(T(pi) ./ T(2) .+ T(pi) ./ T(2) .* padarr ./ pad).^2 .* exp.(im .* angle.(A[low:left-1,:]))
   # pass band
   A[left:right-1,:] .= exp.(im .* angle.(A[left:right-1,:]))
   # right tapering
   A[right:high-1,:] .= cos.(T(pi) ./ T(2) .* padarr ./ pad).^2 .* exp.(im .* angle.(A[right:high-1,:]))
   # right zero cut-off
   A[high:end,:] .= compzero
   return nothing
end
whiten(A::AbstractArray, freqmin::Float64, freqmax::Float64, fs::Float64, N::Int;
    pad::Int=50) = (U = deepcopy(A);
    whiten!(U,freqmin,freqmax,fs,N,pad=pad);
    return U)
"""
   whiten(F, freqmin, freqmax)

Whiten spectrum of FFTData `F` between frequencies `freqmin` and `freqmax`.
Uses real fft to speed up computation.
Returns the whitened (single-sided) fft of the time series.

# Arguments
- `F::FFTData`: FFTData object of fft'd ambient noise data.
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
- `pad::Int`: Number of tapering points outside whitening band.
"""
function whiten!(F::FFTData, freqmin::Float64, freqmax::Float64;pad::Int=50)
    if freqmin < F.freqmin && freqmax > F.freqmax
        @warn "Whitening frequencies ($freqmin, $freqmax Hz) are wider than frequencies
        in FFTData ($(F.freqmin),$(F.freqmax) Hz). Whitening in ($(F.freqmin),$(F.freqmax) Hz) band."
    elseif freqmin < F.freqmin
        @warn "Low whitening frequency $freqmin Hz is lower than minumum frequency
        in FFTData ($(F.freqmin) Hz). Whitening in ($(F.freqmin),$freqmax Hz) band."
    elseif freqmax > F.freqmax
        @warn "High whitening frequency $freqmax Hz is higher than maximum frequency
        in FFTData ($(F.freqmax) Hz). Whitening in ($freqmin,$(F.freqmax) Hz) band."
    end

    N = convert(Int, F.fs * F.cc_len) # number of data points
    freqmin = max(freqmin,F.freqmin) # check for freqmin = 0
    freqmax = min(freqmax,max(F.freqmax,1 / F.cc_len)) # check for freqmax = 0
    whiten!(F.fft, freqmin, freqmax, F.fs, N, pad=pad)
    F.whitened = true
    F.freqmin = freqmin
    F.freqmax = freqmax
    return nothing
end
whiten(F::FFTData, freqmin::Float64, freqmax::Float64;pad::Int=50) =
      (U = deepcopy(F); whiten!(U,freqmin,freqmax,pad=pad);return U)

function whiten!(R::RawData,freqmin::Float64, freqmax::Float64; pad::Int=50)
    if freqmin < R.freqmin && freqmax > R.freqmax
        @warn "Whitening frequencies ($freqmin, $freqmax Hz) are wider than frequencies
        in RawData ($(R.freqmin),$(R.freqmax) Hz). Whitening in ($(R.freqmin),$(R.freqmax) Hz) band."
    elseif freqmin < R.freqmin
        @warn "Low whitening frequency $freqmin Hz is lower than minumum frequency
        in RawData ($(R.freqmin) Hz). Whitening in ($(R.freqmin),$freqmax Hz) band."
    elseif freqmax > R.freqmax
        @warn "High whitening frequency $freqmax Hz is higher than maximum frequency
        in FFTData ($(R.freqmax) Hz). Whitening in ($freqmin,$(R.freqmax) Hz) band."
    end

    N = convert(Int, R.fs * R.cc_len) # number of data points
    freqmin = max(freqmin,R.freqmin) # check for freqmin = 0
    freqmax = min(freqmax,max(R.freqmax,1 / R.cc_len)) # check for freqmax = 0
    FFT = rfft(R.x,1)
    whiten!(FFT,freqmin,freqmax,R.fs, N, pad=pad)
    R.x .= irfft(FFT,N,1)
    R.freqmin = freqmin
    R.freqmax = freqmax
    R.whitened = true
    return nothing
end
whiten(R::RawData,freqmin::Float64, freqmax::Float64; pad::Int=50) =
      (U = deepcopy(R); whiten!(U,freqmin,freqmax,pad=pad);return U)

"""

  coherence!(F,half_win, water_level)

Apply coherence method to FFTData `F`. Where,
``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_A(ω) ∣  ∣ u_B(ω) ∣}``

# Arguments
- `F::FFTData`: FFTData object of fft'd ambient noise data.
- `half_win::Int`: Number of points in half-window to smooth spectrum.
- `water_level::AbstractFloat`: Regularization parameter for spectral smoothing.
                                0.01 is a common value [Mehta, 2007].
"""
function coherence!(F::FFTData, half_win::Int,
                    water_level::Union{Nothing,AbstractFloat}=nothing)
    smoothF = smooth(abs.(F.fft),half_win)
    if !isnothing(water_level)
        reg = water_level .* mean(abs.(F.fft),dims=1)
        smoothF .+= reg
    end
    F.fft ./= smoothF
    return nothing
end
coherence(F::FFTData,half_win::Int,
          water_level::Union{Nothing,AbstractFloat}=nothing) =
          (U = deepcopy(F);coherence!(U,half_win,water_level);
          return U)

"""

  deconvolution!(F,half_win, water_level)

Apply deconvolution method to FFTData `F`. Where,
``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_B(ω) ∣^2}``

# Arguments
- `F::FFTData`: FFTData object of fft'd ambient noise data.
- `half_win::Int`: Number of points in half-window to smooth spectrum.
- `water_level::AbstractFloat`: Regularization parameter for spectral smoothing.
                              0.01 is a common value [Mehta, 2007].
"""
function deconvolution!(F::FFTData, half_win::Int,
                        water_level::Union{Nothing,AbstractFloat}=nothing)
    smoothF = smooth(abs.(F.fft).^2,half_win)
    if !isnothing(water_level)
        reg = water_level .* mean(abs.(F.fft).^2,dims=1)
        smoothF .+= reg
    end
    F.fft ./= smoothF
    return nothing
end
deconvolution(F::FFTData,half_win::Int,
              water_level::Union{Nothing,AbstractFloat}=nothing) =
              (U = deepcopy(F);deconvolution!(U,half_win,water_level);
              return U)

"""
  corrmap(A,maxlag,corr_type,OUTDIR)

Compute cross-correlations using a parallel map.

`corrmap` takes in an array of FFT's and cross-correlates using full parallism.
For a list of `N` `FFT`s, there are `N * (N - 1) / 2` total possible
cross-correlations. To avoid parallel I/O issues, `corrmap` copies the current
correlation `F`, then maps the correlation of `F` against all other correlations
in parallel.

# Arguments
- `A::AbstractArray`: Array of FFTData objects.
- `maxlag::Float64`: Maximum lag time (in seconds) in cross-correlation to save,
                   e.g. `maxlag = 20.` will save lag times = -20.:20. s.
- `corr_type::String`: Type of correlation: `CC` or `PCC`.
- `smoothing_half_win::Int`: Number of points to smooth spectrum of
                  cross-correlations if using deconvolution or coherence.
- `smooth_type::String`: Type of smoothing to apply: `cross-correlation`, `coherence` or
                 `deconvolution`.
- `OUTDIR::String`: Path to save correlation, e.g. "/home/ubuntu/CORR/".
"""
function corrmap(A::Array{FFTData,1},maxlag::Float64,OUTDIR::String;
               corr_type::String="CC",interval::DatePeriod=Day(0),
               smooth_type::String="none",smoothing_half_win::Int=5,
               water_level::Union{Nothing,AbstractFloat}=nothing)
    N = size(A,1)

    if smooth_type == "coherence"
        A = pmap(coherence,A,fill(smoothing_half_win,N),fill(water_level,N))
    end

    # copy the current FFT and correlate against all remaining
    for ii = 1:N-1
        FFT = pop!(A)
        out = fill(deepcopy(FFT),length(A))
        if smooth_type == "deconvolution"
            deconvolution!(FFT,smoothing_half_win,water_level)
        end
        pmap(map_cc,out,A,fill(maxlag,length(A)),
        fill(corr_type,length(A)),fill(OUTDIR,length(A)),
        fill(interval,length(A)))
    end
end


"""
map_cc(FFT1,FFT2,maxlag,OUTDIR)

Input function for corrmap.

Correlates `FFT1` and `FFT2`.

# Arguments
- `FFT1::FFTData`: FFTData object of fft'd ambient noise data.
- `FFT2::FFTData`: FFTData object of fft'd ambient noise data.
- `maxlag::Float64`: Maximum lag time (in seconds) in cross-correlation to save,
               e.g. `maxlag = 20.` will save lag times = -20.:20. s.
- `smooth_type::String`: Type of smoothing to apply: `cross-correlation`, `coherence` or
                 `deconv`.
- `OUTDIR::String`: Path to save correlation, e.g. "/home/ubuntu/CORR/".
-`interval::DatePeriod`:

"""
function map_cc(FFT1::FFTData,FFT2::FFTData,maxlag::Float64,
              corr_type::String,OUTDIR::String,
              interval::DatePeriod)
    println("Correlation $(FFT1.name), $(FFT2.name)")
    C = correlate(FFT1,FFT2,maxlag,corr_type=corr_type)
    if interval > Second(0)
        stack!(C,interval=interval)
    end
    save_corr(C,OUTDIR)
    return nothing
end
