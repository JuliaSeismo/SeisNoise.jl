export process_raw, process_raw!

"""
    compute_fft()


Computes windowed fft of ambient noise data.

Cross-correlates data using either cross-correlation, deconvolution or
cross-coherence. Saves cross-correlations in JLD2 data set.

TO DO:
    - load in data []
    - process_raw [X]
    - check start/end times []
    - chop into matrix []
    - normalize time / freq domain []
    - take fft (with plan_fft) []
    - get parameters for each window []
    - save fft and parameters to JLD2 []


:type maxlag: int
:param maxlag: maximum lag, in seconds, in cross-correlation
:type downsamp_freq: float
:param downsamp_freq: Frequency to which waveforms in stream are downsampled
:return: Downsampled trace or stream object
:type min_dist: float
:param min_dist: minimum distance between stations in km
:type max_dist: float
:param max_dist: maximum distance between stations in km
:type freqmin: float
:param freqmin: minimun frequency for whitening
:type freqmax: float
:param freqmax: maximum frequency for whitening
:type step: float
:param step: time, in seconds, between success cross-correlation windows
:type step: float
:param step: length of noise data window, in seconds, to cross-correlate
"""
function compute_fft()
end

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
            process_raw!(U,fs); return U)

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
      process_raw!(S[i],fs)
  end
  return nothing
end
process_raw(S::SeisData, fs::Real) = (U = deepcopy(S);
            process_raw!(U,fs); return U)

"""
    pre_process(args)

apply 1-bit, filter, whitening
"""
function pre_process()

end


"""

Chop seischannel into matrix with windowing based on noise window and overlap
"""
function window_seis(C::SeisChannel)

end

"""
    remove_resp(args)

remove instrument response - will require reading stationXML and extracting poles
and zeros
"""
function remove_resp(args)
end
