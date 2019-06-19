export phase_shift, phase_shift!

"""
    phase_shift!(C::SeisChannel)

Phase shift SeisChannel if starttime is not aligned with sampling rate.
"""
function phase_shift!(C::SeisChannel)

    # get time offset from sampling rate
    C.t[1,2] = round(C.t[1,2],sigdigits=13)
    t = C.t[1,2]
    dt = 1. / C.fs
    off = mod(millisecond(u2d(t*1e-6)) * 1e-3, dt)
    n = length(C.x)

    if dt - off <= eps(Float64)
        off = 0.
    end

    if off != 0.
        if off <= dt / 2.
            off = -off
        else
            off = dt - off
        end

        # phase shift data
        freq = rfftfreq(length(C.x),C.fs)
        fftdata = rfft(C.x)
        fftdata .= fftdata .* exp.(1im .* 2 .* pi .* freq .* off)
        C.x[:] = irfft(fftdata,length(C.x))

        # reset time of signal
        C.t[1,2] += off * 1e6
    end
    return nothing
end
phase_shift(C::SeisChannel) = (U = deepcopy(C); phase_shift!(U);return U)

"""
    phase_shift!(S::SeisData)

Phase shift SeisData if starttime is not aligned with sampling rate.
"""
function phase_shift!(S::SeisData)
    @inbounds for i = 1:S.n
        phase_shift!(S[i])
    end
    return nothing
end
phase_shift(S::SeisData) = (U = deepcopy(S);phase_shift!(U);return U)
