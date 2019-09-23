export phase_shift, phase_shift!

"""
  phase_shift!(C::SeisChannel)

Phase shift SeisChannel if starttime is not aligned with sampling period.

For example, a channel with starttime of 2016-07-04T00:00:00.008 and a sampling
rate of 100 Hz, the starttime is off the sampling period (0.01) by 2 millisecond
and should be at 2016-07-04T00:00:00.01. The waveform will be phase shifted by
2 milliseconds in the frequency domain. If ϕshift = false, will only correct
starttime of channel.

# Arguments
- `C::SeisChannel`: SeisChannel.
- `ϕshift::Bool`: Boolean to allow/disable phase-shifting.
"""
function phase_shift!(C::SeisChannel; ϕshift::Bool=true)

    # get time offset from sampling rate
    dt = Float32(1. / C.fs)
    off = mod(mod(C.t[1,2],1000000)*Float32(1e-6),dt)
    n = length(C.x)

    if dt - off <= eps(eltype(off))
        off = 0.
    end

    if off != 0.
        if off <= dt / 2.
            off = -off
        else
            off = dt - off
        end

        # phase shift data
        if ϕshift
            freq = rfftfreq(length(C.x),C.fs)
            fftdata = rfft(C.x)
            fftdata .= fftdata .* exp.(1im .* 2 .* pi .* freq .* off)
            C.x[:] = irfft(fftdata,length(C.x))
        end

        # reset time of signal
        C.t[1,2] += off * 1e6
    end
    return nothing
end
phase_shift(C::SeisChannel;ϕshift::Bool=true) = (U = deepcopy(C);
                phase_shift!(U,ϕshift=ϕshift);return U)

"""

  phase_shift!(S::SeisData)

For example, a channel with starttime of 2016-07-04T00:00:00.008 and a sampling
rate of 100 Hz, the starttime is off the sampling period (0.01) by 2 millisecond
and should be at 2016-07-04T00:00:00.01. The waveform will be phase shifted by
2 milliseconds in the frequency domain. If ϕshift = false, will only correct
starttime of channel.

This function will not work properly for sampling rates >= 1000 Hz due to Julia
DateTime resolving time only to milliseconds!

# Arguments
- `S::SeisData`: SeisData.
- `ϕshift::Bool`: Boolean to allow/disable phase-shifting.
"""
function phase_shift!(S::SeisData;ϕshift::Bool=true)
    @inbounds for i = 1:S.n
        phase_shift!(S[i],ϕshift=ϕshift)
    end
    return nothing
end
phase_shift(S::SeisData;ϕshift::Bool=true) = (U = deepcopy(S);
            phase_shift!(U,ϕshift=ϕshift);return U)
