export cwt, icwt, wts_dvv

# """
#
#   wtdtw(ref,cur,t,window,freqmin,freqmax)
#
# dv/v with dynamic time warping method from continuous wavelet transformation.
#
# This function uses the dynamic time warping method at each frequency from the
# continuous wavelet transform to compare the Reference waveform to the current
# waveform to get the relative seismic velocity variation (and associated error).
#
# # Arguments
# - `ref::AbstractArray`: Reference correlation.
# - `cur::AbstractArray`: Current correlation.
# - `t::AbstractArray`: time vector, common to both `ref` and `cur`.
# - `window::AbstractArray`: vector of the indices of the `cur` and `ref` windows
#                           on which you want to do the measurements
# - `freqmin::Float64`: minimum frequency in the correlation [Hz]
# - `freqmax::Float64`: maximum frequency in the correlation [Hz]
# - `maxlag::Int64`: number of maxlag id to search the distance.
# - `b::Int64`: value to control in distance calculation algorithm (see Mikesell et al. 2015).
# - `direction::Int64`: direction to accumulate errors (1=forward, -1=backward, 0=double to smooth)
# - `f0::Real`: Nondimensional frequency from Torrence & Campo, 1998 eq. 1 [Hz].
# - `dj::AbstractFloat`: Spacing between discrete scales. Default value is 1/12.
# - `standardize::Bool`: Remove mean and std from wavelet spectrum or not.
#
# # Returns
# - `dvv::AbstactArray`: Relative Velocity Change dv/v (in %)
# - `err::AbstractArray`: Errors in the dv/v measurements
#
# Originally written in python by Congcong Yuan (30 Jun, 2019)
# """
# function wtdtw(ref,cur,t,window,freqmin,freqmax;maxlag=80,b=1,direction=1,f0=6,dj=1/12,
#                standardize=true)
#
#     T = eltype(ref)
#     dt = mean(diff(t))
#     # apply cwt on two traces
#     cwt1,sj,freqs,coi = cwt(ref,dt,freqmin,freqmax, f0=f0,dj=dj)
#     cwt2,sj,freqs,coi = cwt(cur,dt,freqmin,freqmax, f0=f0,dj=dj)
#
#     # extract real values of cwt
#     rcwt1, rcwt2 = real.(cwt1), real.(cwt2)
#
#     # Use DTW method to extract dvv
#     Nfreq = length(freqs)
#     dvv = zeros(T,Nfreq)
#     err = zeros(T,Nfreq)
#
#     if standardize
#         standardize!(rcwt1)
#         standardize!(rcwt2)
#     end
#
#     for ii = 1:Nfreq
#
#         # find cone of influence
#         indcoi = findall(freqs[ii] .>= 1. ./coi)
#         indt = findall(x->x in indcoi, window)
#
#         if length(indt) == 0
#             continue
#         end
#
#         dvv[ii], err[ii] = dtw_dvv(rcwt1[window[indt],ii], rcwt2[window[indt],ii],t[window[indt]],maxlag, b, direction)
#     end
#     return dvv, err, freqs
# end

"""

  wts_dvv(ref,cur,t,window,freqmin,freqmax)

dv/v with stretching method from continuous wavelet transformation.

This function uses the stretching method at each frequency from the
continuous wavelet transform to compare the Reference waveform to the current
waveform to get the relative seismic velocity variation (and associated error).

# Arguments
- `ref::AbstractArray`: Reference correlation.
- `cur::AbstractArray`: Current correlation.
- `t::AbstractArray`: time vector, common to both `ref` and `cur`.
- `window::AbstractArray`: vector of the indices of the `cur` and `ref` windows
                          on which you want to do the measurements
- `freqmin::Float64`: minimum frequency in the correlation [Hz]
- `freqmax::Float64`: maximum frequency in the correlation [Hz]
- `f0::Real`: Nondimensional frequency from Torrence & Campo, 1998 eq. 1 [Hz].
- `dj::AbstractFloat`: Spacing between discrete scales. Default value is 1/12.
- `standardize::Bool`: Remove mean and std from wavelet spectrum or not.

# Returns
- `dvv::AbstactArray`: Relative Velocity Change dv/v (in %)
- `err::AbstractArray`: Errors in the dv/v measurements

Originally written in python by Congcong Yuan (30 Jun, 2019)
"""
function wts_dvv(ref,cur,t,window,freqmin,freqmax;f0=6,dj=1/12,
               standardize=true)

    T = eltype(ref)
    dt = mean(diff(t))
    # apply cwt on two traces
    cwt1,sj,freqs,coi = cwt(ref,dt,freqmin,freqmax, f0=f0,dj=dj)
    cwt2,sj,freqs,coi = cwt(cur,dt,freqmin,freqmax, f0=f0,dj=dj)

    # Use DTW method to extract dvv
    Nfreq = length(freqs)
    dvv = zeros(T,Nfreq)
    err = zeros(T,Nfreq)

    for ii = 1:Nfreq
        icwt1 = icwt(cwt1[:,ii],sj[ii],dt)
        icwt2 = icwt(cwt2[:,ii],sj[ii],dt)

        if standardize
            standardize!(icwt1)
            standardize!(icwt2)
        end

        # find cone of influence
        indcoi = findall(freqs[ii] .>= 1. ./coi)
        indt = findall(x->x in indcoi, window)

        if length(indt) == 0
            continue
        end

        # calculate dvv with stretching
        dvv[ii], cc,cdp, ep, err[ii],allC = stretching(icwt1, icwt2,t,window[indt],
               freqmin,freqmax,ntrial=30,dvmin=-0.01,dvmax=0.01)
    end
    return dvv, err, freqs
end


# """
#
#   dtw_dvv(ref,cur,t,maxlag,d,direction)
#
# dv/v with dynamic time warping method.
#
# This function uses the dynamic time warping method to compare the optimal stretching
# between Reference waveform to the current waveform to get the relative seismic
# velocity variation (and associated error).
#
# # Arguments
# - `ref::AbstractArray`: Reference correlation.
# - `cur::AbstractArray`: Current correlation.
# - `t::AbstractArray`: time vector, common to both `ref` and `cur`.
# - `maxlag::Int64`: number of maxlag id to search the distance.
# - `b::Int64`: value to control in distance calculation algorithm (see Mikesell et al. 2015).
# - `direction::Int64`: direction to accumulate errors (1=forward, -1=backward, 0=double to smooth)
#
# # Returns
# - `dvv::Float64`: Relative Velocity Change dv/v (in %)
# - `err::Float64`: Errors in the dv/v measurements
#
# """
# function dtw_dvv(ref,cur,t,maxlag, b, direction)
#     dt = mean(diff(t))
#     stbarTime, stbar, dist, error = dtwdt(ref, cur, dt, maxLag=maxlag, b=b, direction=direction)
#
#     # perform linear regression
#     model = glm(@formula(Y ~0 + X),DataFrame(X=t,Y=stbarTime),Normal(),
#                 IdentityLink(),wts=ones(length(t)))
#
#     return coef(model)[1], stderror(model)[1]
# end

"""

  cwt(signal,dt,freqmin,freqmax)

Continuous wavelet transform of the signal at specified scales.
Note: uses the Morlet wavelet ONLY.

# Arguments
- `signal::AbstractArray`: Input signal array.
- `dt::AbstractFloat`: Sampling interval [s].
- `freqmin::AbstractFloat`: Minimum frequency for cwt [Hz].
- `freqmax::AbstractFloat`: Maximum frequency for cwt [Hz].
- `f0::Real`: Nondimensional frequency from Torrence & Campo, 1998 eq. 1 [Hz].
- `dj::AbstractFloat`: Spacing between discrete scales. Default value is 1/12.

# Returns
- `W::AbstractArray`: Wavelet transform from Morlet wavelet.
- `freqs::AbstractArray`: Fourier frequencies for wavelet scales [Hz].
- `coi::AbstractArray`: Cone of influence - maximum period (in s) of useful
    information at that particular time. Periods greater than
    those are subject to edge effects.

This is a Julia translation of the cwt in pycwt https://github.com/regeirk/pycwt
"""
function cwt(signal::AbstractArray{T,1},dt::AbstractFloat,freqmin::AbstractFloat,
             freqmax::AbstractFloat;f0=6.,dj=1/12) where T <: AbstractFloat
    n0 = length(signal)
    flambda = T.((4 .* π) ./ (f0 .+ sqrt(2 .+ f0^2)))
    s0 = 2 * dt / flambda
    J = convert(Int,round(log2(n0 .* dt ./ s0) ./ dj))

    # The scales as of Mallat 1999
    sj = T.(s0 .* 2 .^ ((0:J) .* dj))
    freqs = 1 ./ (flambda .* sj)

    # subset by freqmin & freqmax
    ind = findall((freqs .> freqmin) .& (freqs .< freqmax))
    sj = sj[ind]
    freqs = freqs[ind]

    # signal fft
    signal_ft = fft(signal,1)
    N = length(signal_ft)
    # Fourier angular frequencies
    ftfreqs = T.(2 .* π * FFTW.fftfreq(n0,1/dt))

    # Creates wavelet transform matrix as outer product of scaled transformed
    # wavelets and transformed signal according to the convolution theorem.
    # (i)   Transform scales to column vector for outer product;
    # (ii)  Calculate 2D matrix [s, f] for each scale s and Fourier angular
    #       frequency f;
    # (iii) Calculate wavelet transform;
    psi_ft_bar = ((sj' .* ftfreqs[2] .* T(N)) .^ T(0.5)) .* conj.(psi_ft(ftfreqs * sj',f0))
    W = ifft(signal_ft .* psi_ft_bar,1)

    # Checks for NaN in transform results and removes them from the scales if
    # needed, frequencies and wavelet transform. Trims wavelet transform at
    # length `n0`.
    sel = findall(.!all(isnan.(W),dims=1)[:])
    if length(sel) > 0
        sj = sj[sel]
        freqs = freqs[sel]
        W = W[:,sel]
    end

    # Determines the cone-of-influence. Note that it is returned as a function
    # of time in Fourier periods. Uses triangualr Bartlett window with
    # non-zero end-points.
    coi = n0 ./ 2 .- abs.((0:n0-1) .- (n0 - 1.) ./2.)
    coi = T.(flambda ./ sqrt(2.) .* dt .* coi)
    return W, sj, freqs, coi
end

"""

  icwt(W,sj,dt)

Inverse continuous wavelet transform at specified scales.
Note: uses the Morlet wavelet ONLY.

# Arguments
- `W::AbstractArray`: Wavelet transform, the result of the `cwt` function.
- `sj::AbstractArray`: Scale indices as returned by the `cwt` function.
- `dt::AbstractFloat`: Sampling interval [s].
- `dj::AbstractFloat`: Spacing between discrete scales. Default value is 1/12.

# Returns
- `iW::AbstractArray`: Inverse wavelet transform from Morlet wavelet.

This is a Julia translation of the icwt in pycwt https://github.com/regeirk/pycwt
Note that the pycwt version has incorrect scaling.
"""
function icwt(W::AbstractArray, sj::AbstractArray, dt::AbstractFloat;dj=1/12)
    T = real(eltype(W))
    # As of Torrence and Compo (1998), eq. (11)
    iW = T(dj .* sqrt(dt) ./ 0.776 .* (pi ^ 0.25)) .* sum(real.(W) ./ (sj .^ T(0.5))',dims=2)
    return iW
end

"""

  icwt(W,sj,dt)

Inverse continuous wavelet transform at specified scale.
Note: uses the Morlet wavelet ONLY.

# Arguments
- `W::AbstractArray`: Wavelet transform, the result of the `cwt` function.
- `sj::AbstractFloat`: Scale index as returned by the `cwt` function.
- `dt::AbstractFloat`: Sampling interval [s].
- `dj::AbstractFloat`: Spacing between discrete scales. Default value is 1/12.

# Returns
- `iW::AbstractArray`: Inverse wavelet transform from Morlet wavelet.

This is a Julia translation of the icwt in pycwt https://github.com/regeirk/pycwt
Note that the pycwt version has incorrect scaling.
"""
function icwt(W::AbstractArray, sj::AbstractFloat, dt::AbstractFloat;dj=1/12)
    T = real(eltype(W))
    # As of Torrence and Compo (1998), eq. (11)
    iW = T(dj .* sqrt(dt) ./ 0.776 .* (pi ^ 0.25)) .* real.(W) ./ sj ^ T(0.5)
    return iW
end

function psi_ft(A::AbstractArray{T},f0::Real) where T <: AbstractFloat
    return exp.(T(-0.5) .* (A .- T(f0)) .^2) .* T(π ^ -0.25)
end
