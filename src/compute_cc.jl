export process_raw, process_raw!

"""
    process_raw!(C,fs)

Pre-process month-long stream of data.

Checks:
- sample rate is fs
- downsamples data
- checks for gaps in data
- phase-shifts data to begin at 00:00:00.0

# Arguments
- `C::SeisChannel`: SeisChannel structure.
- `fs::Real`: Sampling rate to downsample `S`.
"""
function process_raw!(C::SeisChannel, fs::Real)
    demean!(C)        # remove mean from channel
    ungap!(C)         # replace gaps with mean of channel
    detrend!(C)       # remove linear trend from channel
    taper!(C)         # taper channel ends
    lowpass!(C,fs)    # lowpass filter before downsampling
    C = downsample!(C,fs) # downsample to lower fs
    check_and_phase_shift!(C) # timing offset from sampling period
    return nothing
end
process_raw(C::SeisChannel, fs::Real) = (U = deepcopy(C);
            process_raw!(U); return U)

"""
  process_raw!(S::SeisData, fs::Real)

  Pre-process month-long stream of data.

  Checks:
  - sample rate is fs
  - downsamples data
  - checks for gaps in data
  - phase-shifts data to begin at 00:00:00.0

  # Arguments
  - `C::SeisChannel`: SeisChannel structure.
  - `fs::Real`: Sampling rate to downsample `S`.
"""
function process_raw!(S::SeisData, fs::Real)
  @inbounds for i = 1:S.n
      process_raw!(S[i])
  end
  return nothing
end
process_raw(S::SeisData, fs::Real) = (U = deepcopy(S);
            process_raw!(U,fs); return U)

"""
    check_and_phase_shift!(C::SeisChannel)

Phase shift SeisChannel if starttime is not aligned with sampling rate.
"""
function check_and_phase_shift!(C::SeisChannel)
    t = C.t[1,2]
    dt = 1. / C.fs
    off = mod(millisecond(u2d(t)) * 1e-3, 1. / C.fs)
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
        nfft = nextprod([2, 3, 5],n)
        C.x[:] = [C.x; zeros(eltype(C.x), nfft - n)]
        freq = fftfreq(nfft,C.fs)
        fftdata = fft(C.x)
        fftdata .= fftdata .* exp.(1im .* 2 .* pi .* fftfreq .* dt)
        C.x[:] = ifft(fftdata)[1:n]
        C.t[1,2] += off * 1e6
    end
    return nothing
end
check_and_phase_shift(C::SeisChannel) = (U = deepcopy(C);
                                         check_and_phase_shift!(U);
                                         return U)

"""
    check_and_phase_shift!(S::SeisData)

Phase shift SeisData if starttime is not aligned with sampling rate.
"""
function check_and_phase_shift!(S::SeisData)
    @inbounds for i = 1:S.n
        check_and_phase_shift!(S[i])
    end
    return nothing
end
check_and_phase_shift(S::SeisData) = (U = deepcopy(S);
                                      check_and_phase_shift!(U);
                                      return U)
