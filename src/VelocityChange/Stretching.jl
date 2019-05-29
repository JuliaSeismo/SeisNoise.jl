export stretching

"""

  stretching(ref,cur,t,window,fmin,fmax;dvmin,dvmax,ntrials)

This function compares the Reference waveform to stretched/compressed current
waveforms to get the relative seismic velocity variation (and associated error).
It also computes the correlation coefficient between the Reference waveform and
the current waveform.

# Arguments
- `ref::AbstractArray`: Reference correlation.
- `cur::AbstractArray`: Current correlation.
- `t::AbstractArray`: time vector, common to both `ref` and `cur`.
- `window:AbstractArray`: vector of the indices of the `cur` and `ref` windows
                          on which you want to do the measurements
- `fmin:Float64`: minimum frequency in the correlation [Hz]
- `fmax:Float64`: maximum frequency in the correlation [Hz]
- `dvmin:Float64`: minimum bound for the velocity variation; e.g. dvmin=-0.03
                   for -3% of relative velocity change
- `dvmax:Float64`: maximum bound for the velocity variation; e.g. dvmin=0.03
                  for 3% of relative velocity change
- `ntrial`:  number of stretching coefficient between dvmin and dvmax, no need to be higher than 100

# Returns
- `dv:AFloat64`: Relative Velocity Change dv/v (in %)
- `cc:Float64`: Correlation coefficient between the reference waveform and the
                      best stretched/compressed current waveform
- `cdp:Float64`: Correlation coefficient between the reference waveform and the
                 initial current waveform
- `ϵ:Array{Float64,1}`: Vector of Epsilon values (ϵ =-dt/t = dv/v)
- `err:Float64`: Errors in the dv/v measurements based on [Weaver et al., 2011](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-246X.2011.05015.x)
- `C:Array{Float64,1}`: Vector of the correlation coefficient between the
                        reference waveform and every stretched/compressed
                        current waveforms

This code is a Julia translation of the Python code from [Viens et al., 2018](https://github.com/lviens/2018_JGR).
"""
function stretching(ref::AbstractArray,cur::AbstractArray,t::AbstractArray,
                    window::AbstractArray,fmin::Float64,fmax::Float64;
                    dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial::Int=100)
    ϵ = range(dvmin, stop=dvmax, length=ntrial)
    L = 1. .+ ϵ
    tau = t * L'
    C = zeros(ntrial)

    # set of stretched/compressed current waveforms
    waveform_ref = ref[window]
    for ii = 1:ntrial
        s = LinearInterpolation(tau[:,ii],cur,extrapolation_bc=Flat())(t)
        waveform_ref = ref[window]
        waveform_cur = s[window]
        C[ii] = cor(waveform_ref,waveform_cur)
    end

    cdp = cor(cur[window],ref[window])
    # find the maximum correlation coefficient
    imax = argmax(C)
    if imax >= ntrial -1
        imax = ntrial - 2
    end
    if imax <= 3
        imax += 3
    end
    # Proceed to the second step to get a more precise dv/v measurement
    dtfiner = Array(range(ϵ[imax-2],stop=ϵ[imax+1],length=500))
    # to get same result as scipy, use line below; this requires gcc
    # using Dierckx; etp = Spline1D(ϵ[range(imax-3,stop=imax+1)],C[range(imax-3,stop=imax+1)])
    etp = CubicSplineInterpolation(ϵ[range(imax-3,stop=imax+1)],C[range(imax-3,stop=imax+1)])
    CCfiner = etp(dtfiner)
    dv = 100. * dtfiner[argmax(CCfiner)]
    cc = maximum(CCfiner) # Maximum correlation coefficient of the refined analysis

    # Error computation based on Weaver, R., C. Hadziioannou, E. Larose, and M.
    # Campillo (2011), On the precision of noise-correlation interferometry,
    # Geophys. J. Int., 185(3), 1384?1392
    T = 1 / (fmax - fmin)
    X = cc
    wc = π * (fmin + fmax)
    tmin = t[window][1]
    tmax = t[window][end]
    t1 = minimum([tmin,tmax])
    t2 = maximum([tmin,tmax])
    err = 100 * (sqrt(1-X^2)/(2*X)*sqrt((6*sqrt(π/2)*T)/(wc^2*(t2^3-t1^3))))

    return dv,cc,cdp,Array(ϵ),err,C
end
