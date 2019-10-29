# cross-correlation module
export clean_up!, clean_up, correlate, compute_cc, generate_pairs, map_cc, corrmap
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
                   corners::Int=4, zerophase::Bool=true)
    demean!(A)
    detrend!(A)
    taper!(A,fs)
    bandpass!(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
    return nothing
end
clean_up(A::AbstractArray, freqmin::Real, freqmax::Real, fs::Real;
         corners::Int=4, zerophase::Bool=false) =
                 (U = deepcopy(A); clean_up!(U,freqmin,freqmax, fs,
                  corners=corners, zerophase=zerophase); return U)

clean_up!(C::CorrData,freqmin::Float64,freqmax::Float64; corners::Int=4,
          zerophase::Bool=true) = clean_up!(C.corr,freqmin,freqmax,C.fs,
          corners=corners,zerophase=zerophase)

clean_up!(R::RawData,freqmin::Float64,freqmax::Float64; corners::Int=4,
          zerophase::Bool=true) = clean_up!(C.x,freqmin,freqmax,R.fs,
          corners=corners,zerophase=true)

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
function correlate(FFT1::AbstractArray, FFT2::AbstractArray, N::Int, maxlag::Int)

    corrF = conj.(FFT1) .* FFT2

    # take inverse fft
    corrT = irfft(corrF,N,1)
    corrT = fftshift(corrT,1)

    # return corr[-maxlag:maxlag]
    t = range(-Int(N/2) + 1, stop=Int(N/2) - 1)
    ind = findall(abs.(t) .<= maxlag)
    corrT = corrT[ind,:]
end

"""
    compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;
               corr_type::String="cross-correlation")

Cross-correlate ambient noise data in the frequency domain.

Cross-correlation can be done using one of three options:

- Cross-correlation: ``C_{AB}(ω) = u_A(ω) u^∗_B(ω)``
- Coherence: ``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_A(ω) ∣  ∣ u_B(ω) ∣}``
- Deconvolution: ``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_B(ω) ∣^2}``

# Arguments
- `FFT1::FFTData`: FFTData object of fft'd ambient noise data.
- `FFT2::FFTData`: FFTData object of fft'd ambient noise data.
- `maxlag::Float64`: Maximum lag time (in seconds) in cross-correlation to save,
                     e.g. `maxlag = 20.` will save lag times = -20.:20. s.
- `corr_type::String`: Type of correlation: `cross-correlation`, `coherence` or
                       `deconv`.
"""
function compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;
                    corr_type::String="cross-correlation")

    N = convert(Int,round(FFT1.cc_len * FFT1.fs)) # number of data points
    comp = FFT1.name[end] * FFT2.name[end]
    # get intersect of dates; return nothing if no intersect
    inter = intersect(FFT1.t,FFT2.t)
    if length(inter) == 0
        return nothing
    end

    ind1 = findall(x -> x ∈ inter, FFT1.t)
    ind2 = findall(x -> x ∈ inter, FFT2.t)

    corr = correlate(FFT1.fft[:,ind1], FFT2.fft[:,ind2], N,
                     convert(Int,round(maxlag * FFT1.fs)))
    rotated = false

    return CorrData(FFT1, FFT2, comp, rotated, corr_type,
                    maxlag, inter, corr)

end

"""
    corrmap(A,maxlag,smoothing_half_win,corr_type,OUTDIR)

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
- `smoothing_half_win::Int`: Number of points to smooth spectrum of
                    cross-correlations if using deconvolution or coherence.
- `corr_type::String`: Type of correlation: `cross-correlation`, `coherence` or
                       `deconv`.
- `OUTDIR::String`: Path to save correlation, e.g. "/home/ubuntu/CORR/".
"""
function corrmap(A::Array{FFTData,1},maxlag::Float64, smoothing_half_win::Int,
                 corr_type::String,OUTDIR::String;
                 water_level::Union{Nothing,AbstractFloat}=nothing)
        N = size(A,1)

        if corr_type == "coherence"
            A = pmap(coherence,A,fill(smoothing_half_win,N),fill(water_level,N))
        end

        # copy the current FFT and correlate against all remaining
        for ii = 1:N-1
            FFT = pop!(A)
            if corr_type == "deconv"
                deconvolution!(FFT,smoothing_half_win,water_level)
            end
            out = fill(deepcopy(FFT),length(A))
            pmap(map_cc,out,A,fill(maxlag,length(A)),
                 fill(corr_type,length(A)),fill(OUTDIR,length(A)))
        end
end

function corrmap(A::Array{FFTData,1},maxlag::Float64, smoothing_half_win::Int,
                 corr_type::String,OUTDIR::String,
                 interval::Union{Month,Day,Hour,Second};
                 water_level::Union{Nothing,AbstractFloat}=nothing)
        N = size(A,1)

        if corr_type == "coherence"
            A = pmap(coherence,A,fill(smoothing_half_win,N),fill(water_level,N))
        end

        # copy the current FFT and correlate against all remaining
        for ii = 1:N-1
            FFT = pop!(A)
            out = fill(deepcopy(FFT),length(A))
            if corr_type == "deconv"
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
- `corr_type::String`: Type of correlation: `cross-correlation`, `coherence` or
                   `deconv`.
- `OUTDIR::String`: Path to save correlation, e.g. "/home/ubuntu/CORR/".

"""
function map_cc(FFT1::FFTData,FFT2::FFTData,maxlag::Float64,
                corr_type::String,OUTDIR::String)
    println("Correlation $(FFT1.name), $(FFT2.name)")
    C = compute_cc(FFT1,FFT2,maxlag,corr_type=corr_type)
    if !isnothing(C)
        save_corr(C,OUTDIR)
    end
    return nothing
end

function map_cc(FFT1::FFTData,FFT2::FFTData,maxlag::Float64,
                corr_type::String,OUTDIR::String,
                interval::Union{Month,Day,Hour,Second})
    println("Correlation $(FFT1.name), $(FFT2.name)")
    C = compute_cc(FFT1,FFT2,maxlag,corr_type=corr_type)
    if !isnothing(C)
        stack!(C,interval=interval)
        save_corr(C,OUTDIR)
    end
    return nothing
end

"""
    whiten(A, freqmin, freqmax, fs, pad=50)

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
function whiten!(A::AbstractArray, freqmin::Float64, freqmax::Float64, fs::Float64,
                 N::Int;pad::Int=50)
    T = real(eltype(A))

    # get size and convert to Float32
    Nrows,Ncols = size(A)

    # get whitening frequencies
    freqvec = rfftfreq(N,fs)
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

    # left zero cut-off
    compzero = complex(T(0))
    for jj = 1:Ncols
         for ii = 1:low
            A[ii,jj] = compzero
        end

    # left tapering
         for ii = low+1:left
            A[ii,jj] = cos(T(pi) ./ T(2) .+ T(pi) ./ T(2) .* (ii-low-1) ./ pad).^2 .* exp(
            im .* angle(A[ii,jj]))
        end

    # pass band
         for ii = left+1:right
            A[ii,jj] = exp(im .* angle(A[ii,jj]))
        end

    # right tapering
         for ii = right+1:high
            A[ii,jj] = cos(T(pi) ./ T(2) .* (ii-right) ./ pad).^2 .* exp(
            im .* angle(A[ii,jj]))
        end

    # right zero cut-off
         for ii = high+1:size(A,1)
            A[ii,jj] = compzero
        end
    end
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

    freqmin = max(freqmin,F.freqmin) # check for freqmin = 0
    freqmax = min(freqmax,max(F.freqmax,1 / F.cc_len)) # check for freqmax = 0
    whiten!(F.fft, freqmin, freqmax, F.fs,convert(Int, F.fs * F.cc_len), pad=pad)
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

    freqmin = max(freqmin,R.freqmin) # check for freqmin = 0
    freqmax = min(freqmax,max(R.freqmax,1 / R.cc_len)) # check for freqmax = 0
    FFT = rfft(R.x,1)
    whiten!(FFT,freqmin,freqmax,R.fs, convert(Int, R.fs * R.cc_len), pad=pad)
    R.x .= irfft(FFT,convert(Int,R.cc_len * R.fs),1)
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
end
deconvolution(F::FFTData,half_win::Int,
              water_level::Union{Nothing,AbstractFloat}=nothing) =
              (U = deepcopy(F);deconvolution!(U,half_win,water_level);
              return U)

function generate_pairs(files::AbstractArray)
    N = length(files)
    num_pairs = convert(Int,round(N * (N-1) / 2 + N))
    pairs = Array{Array{String,1},1}(undef,num_pairs)
    count = 0
    for ii = 1:length(files)
      for jj = ii:length(files)
          count += 1
          pairs[count] = [files[ii], files[jj]]
      end
    end
    return pairs
end
