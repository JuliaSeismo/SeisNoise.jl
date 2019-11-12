export mwcs, mwcs_dvv

function getCoherence(dcs,ds1,ds2)
  coh = zeros(size(dcs))
  valids = intersect(findall(x -> abs.(x) > 0,ds1),findall(x -> abs.(x) > 0,ds2))
  coh[valids] = dcs[valids] ./ (ds1[valids] .* ds2[valids])
  coh[coh .> 1.0 ] .= 1.0
  return coh
end

"""

  mwcs(ref, cur, fmin, fmax, fs, tmin, window_length, window_step,
       smoothing_half_win)

Change in velocity measurement using the Moving Window Cross-Spectrum technique.

The current correlation `cur` is compared to the reference correlation `ref`.
Both time series are sliced in several overlapping windows.
Each slice is mean-adjusted and cosine-tapered (85% taper) before being Fourier-
transformed to the frequency domain.
``F_{cur}(ν)`` and ``F_{ref}(ν)`` are the first halves of the
Hermitian symmetric Fourier-transformed segments. The cross-spectrum
``X(ν)`` is defined as
``X(ν) = F_{ref}(ν) F_{cur}^*(ν)``
in which ``{}^*`` denotes the complex conjugation.
``X(ν)`` is then smoothed by convolution with a Hanning window.
The similarity of the two time-series is assessed using the cross-coherency
between energy densities in the frequency domain:
``C(ν) = \\frac{|\\overline{X(ν))}|}{\\sqrt{|\\overline{F_{ref}(ν)|^2} |\\overline{F_{cur}(ν)|^2}}}``
in which the over-line here represents the smoothing of the energy spectra for
``F_{ref}`` and ``F_{cur}`` and of the spectrum of ``X``. The mean
coherence for the segment is defined as the mean of ``C(ν)`` in the
frequency range of interest. The time-delay between the two cross correlations
is found in the unwrapped phase, ``ϕ(ν)``, of the cross spectrum and is
linearly proportional to frequency:
``ϕ_j = m. ν_j, m = 2 π δ t``
The time shift for each window between two signals is the slope ``m`` of a
weighted linear regression of the samples within the frequency band of interest.
The weights are those introduced by [Clarke2011]_,
which incorporate both the cross-spectral amplitude and cross-coherence, unlike
[Poupinet1984]_. The errors are estimated using the weights (thus the
coherence) and the squared misfit to the modelled slope:
``e_m = \\sqrt{\\sum_j{(\\frac{w_j ν_j}{\\sum_i{w_i ν_i^2}})^2}σ_{ϕ}^2}``
where ``w`` are weights, ``ν`` are cross-coherences and
``σ_{ϕ}^2`` is the squared misfit of the data to the modelled slope
and is calculated as ``σ_{ϕ}^2 = \\frac{\\sum_j(ϕ_j - m ν_j)^2}{N-1}``
The output of this process is a table containing, for each moving window: the
central time lag, the measured delay, its error and the mean coherence of the
segment.

# Arguments
- `ref::AbstractArray`: Reference correlation.
- `cur::AbstractArray`: Current correlation.
- `fmin:Float64`: minimum frequency in the correlation [Hz]
- `fmax:Float64`: maximum frequency in the correlation [Hz]
- `fs:Float64`: Sampling rate of `ref` and `cur` [Hz]
- `tmin:Float64`: The leftmost time lag [s]
- `window_len:Float64`: The moving window length [s]
- `window_step:Float64`: The step to jump for the moving window [s]
- `smoothing_half_win:Int`: Defines the half length of the smoothing hanning
                            window.

# Returns
- `time_axis:Array{Float64,1}`: Central time of each moving window [s]
- `dt:Array{Float64,1}`: dt for each moving window
- `err:Array{Float64,1}`: Errors for each moving window
- `mcoh:Array{Float64,1}`: Mean coherence for each moving window

This code is a Julia translation of the Python code from [MSNoise](https://github.com/ROBelgium/MSNoise/blob/master/msnoise/move2obspy.py).
"""
function mwcs(ref::AbstractArray,cur::AbstractArray,fmin::Float64,
                   fmax::Float64,fs::Float64,tmin::Float64,
                   window_length::Float64,window_step::Float64,
                   smoothing_half_win::Int)

    # create time axis for mwcs
    time_axis = Array(tmin + window_length / 2. : window_step : tmin +
                      length(ref) / fs - window_length / 2.)

    window_length_samples = convert(Int,window_length * fs)
    window_step_samples = convert(Int,window_step * fs)
    minind = 1:window_step_samples:length(ref) - window_length_samples
    padd = convert(Int,2 ^ (ceil(log2(abs(window_length_samples))) + 2))

    N = length(minind)
    dt = zeros(N)
    err = zeros(N)
    coh = zeros(N)
    time_axis = time_axis[1:N]

    # Find values in frequency range of interest
    freq_vec = FFTW.rfftfreq(padd,fs)
    index_range = findall(x -> x >= fmin && x <= fmax,freq_vec)
    cci = zeros(window_length_samples,N)
    cri = zeros(window_length_samples,N)

    # fill matrices
    for ii = 1:N
        cci[:,ii] = cur[minind[ii]:minind[ii]+window_length_samples-1]
        cri[:,ii] = ref[minind[ii]:minind[ii]+window_length_samples-1]
    end

    # preprocess
    SeisNoise.demean!(cci)
    SeisNoise.detrend!(cci)
    SeisNoise.taper!(cci,fs,max_percentage=0.425)
    SeisNoise.demean!(cri)
    SeisNoise.detrend!(cri)
    SeisNoise.taper!(cri,fs,max_percentage=0.425)

    # take fourier transform with padding
    fcur = rfft([cci;zeros(eltype(cci),padd-size(cci)[1],N)],1)
    fref = rfft([cri;zeros(eltype(cri),padd-size(cri)[1],N)],1)

    fcur2 = real(fcur).^2 + imag(fcur).^2
    fref2 = real(fref).^2 + imag(fref).^2

    # calculate cross-spectrum
    X = fref .* conj.(fcur)

    if smoothing_half_win != 0
        dcur = sqrt.(smooth(fcur2,smoothing_half_win))
        dref = sqrt.(smooth(fref2,smoothing_half_win))
        X = smooth(X, smoothing_half_win)
    else
        dcur = sqrt.(fcur2)
        dref = sqrt.(fref2)
    end

    dcs = abs.(X)

    # Get coherence and its mean value
    coh = getCoherence(dcs,dref,dcur)
    mcoh = reshape(mean(coh[index_range,:],dims=1),N)

    # Get Weights
    w = 1. ./ (1. ./ coh[index_range,:].^2  .- 1.)
    w[coh[index_range,:] .>= 0.99] .= 1. / (1. / 0.9801 - 1.)
    w = sqrt.(w .* sqrt.(dcs[index_range,:]))
    w = real(w)

    # frequency array
    v = real(freq_vec[index_range]) .* 2π

    # phase
    phi = angle.(X)
    phi[1,:] .= 0.
    unwrap!(phi,dims=1)
    phi = phi[index_range,:]

    # Calculate the slope with a weighted least square linear regression
    # forced through the origin
    for ii = 1:N
        model = glm(@formula(Y ~0 + X),DataFrame(X=v,Y=phi[:,ii]),Normal(),
                    IdentityLink(),wts=w[:,ii])
        # a much faster, unweighted version of this is: m = v \ phi
        dt[ii] = coef(model)[1]

        # Errors
        e = sum((phi[:,ii] .- dt[ii] .* v).^2) / (length(v) - 1)
        s2x2 = sum(v.^2 .* w[:,ii].^2)
        sx2 = sum(w[:,ii] .* v.^2)
        err[ii] = sqrt(e * s2x2 / sx2^2)
    end
    dt ./= 2π
    err ./= 2π
    return time_axis, dt, err, mcoh
end

"""
    mwcs_dvv(time_axis,dt,err,coh,dtt_lag,dist,dtt_v,dtt_minlag,dtt_width,
             dtt_sides,min_coh, max_err, max_dt)

Regresses dv/v from dt/t measurements.

Solves dt = a + bt, where b = dv/v, a = instrumental drift, dt = time lag
at time t. Solves with a weighted linear regression with weights equal to
1/error**2.

# Arguments
- `time_axis:Array{Float64,1}`: Central time of each moving window [s]
- `dt:Array{Float64,1}`: dt for each moving window
- `err:Array{Float64,1}`: Errors for each moving window
- `coh:Array{Float64,1}`: Mean coherence for each moving window
- `dist:Float64`: Distance between stations [km]
- `dtt_lag:String`: Type of time lag 'dynamic' or 'static'. When `dtt_lag`
    is set to "dynamic", the inter-station distance is used to determine
    the minimum time lag
- `dtt_v:Float64`: Velocity for minumum time lag. The velocity is determined by
    the user so that the minlag doesn't include the ballistic waves.
    If ballistic waves are visible with a velocity of 2 km/s, one could
    configure dtt_v=1.0. If stations are located 15 km apart, the minimum
    lag time will be set to 15 s.
- `minlag:Float64`: Statically set minimum lag [s]
- `dtt_width:Float64`: Width of the lag window used [s] i.e., if
    dtt_width = 30, and min_lag = 15, window = 15s - 45s
- `dtt_sides:String`: Sides of cross-correlation function to use. Either
    'Both' or 'left'.
- `min_coh:Float64`: Minimum allowed coherency between reference and current
                     correlation
- `max_err:Float64`: Maximum allowed error in dt/t regression
-  `max_dt:Float64`: Maximum allowed dt from MWCS [s]

# Returns
- `m::Float64`: dt/t for current correlation
- `em::Float64`: Error for calulatoin of `m`
- `a::Float64`: Intercept for regression calculation
- `ea::Float64`: Error on intercept
- `m0::Float64`: dt/t for current correlation with no intercept
- `em0::Float64`: Error for calulatoin of `m0`
"""
function mwcs_dvv(time_axis::AbstractArray, dt::AbstractArray,
                     err::AbstractArray, coh::AbstractArray, dtt_lag::String,
                     dist::Float64, dtt_v::Float64, dtt_minlag::Float64,
                     dtt_width::Float64, dtt_sides::String;
                     min_coh::Float64=0.7, max_err::Float64=1.,
                     max_dt::Float64=1.)

    # lag bins
    if dtt_lag == "static"
        lmlag = -dtt_minlag
        rmlag = dtt_minlag
    else
        lmlag = - dist / dtt_v
        rmlag = dist / dtt_v
    end

    lMlag = lmlag - dtt_width
    rMlag = rmlag + dtt_width

    # subset correlation by sides
    if dtt_sides == "both"
        tindex = findall(x -> abs.(x) .>= rmlag && abs.(x) .<= rMlag, time_axis)
    elseif dtt_sides == "left"
        tindex = findall(x -> x .>= lMlag && x .<= lmlag, time_axis)
    else
        tindex = findall(x -> x .>= rmlag && x .<= rMlag, time_axis)
    end

    # index to subset data
    dt, err, coh = dt[tindex], err[tindex], coh[tindex]
    index = intersect(findall(x -> x .<= min_coh,coh),
                      findall(x -> x .>= max_err,err),
                      findall(x -> abs.(x) .>= max_dt,time_axis))
    dt[index] .= 0.
    err[index] .= 1.
    coh[index] .= 1.

    # weight statistics
    w = 1. ./ err .^ 2

    # regress data using least squares
    VecXfilt = time_axis[tindex]
    model0 = glm(@formula(Y ~0 + X),DataFrame(X=VecXfilt,Y=dt),Normal(),
                IdentityLink(),wts=w)
    model = glm(@formula(Y ~ X),DataFrame(X=VecXfilt,Y=dt),Normal(),
                IdentityLink(),wts=w)

    a,m = coef(model)
    ea, em = stderror(model)

    m0 = coef(model0)[1]
    em0 = stderror(model0)[1]
    return m, em, a, ea, m0, em0
end
