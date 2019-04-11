export process_raw, process_raw!, process_fft, compute_fft

"""
    compute_fft()


Computes windowed fft of ambient noise data.

Cross-correlates data using either cross-correlation, deconvolution or
cross-coherence. Saves cross-correlations in JLD2 data set.

TO DO:
    - load in data []
    - process_raw [x]
    - check start/end times [x]
    - chop into matrix [x]
    - normalize time / freq domain [x]
    - take fft [x]
    - get parameters for each window (amplitude, mad) []
    - save fft and parameters to JLD2 []


:type maxlag: int
:param maxlag: maximum lag, in seconds, in cross-correlation
:type fs: Real
:param fs: Frequency to which waveforms in stream are downsampled
:type freqmin: float
:param freqmin: minimun frequency for whitening
:type freqmax: float
:param freqmax: maximum frequency for whitening
:type cc_step: Real
:param cc_step: time, in seconds, between success cross-correlation windows
:type cc_len: Real
:param cc_len: length of noise data window, in seconds, to cross-correlate
"""
function compute_fft(S::SeisData,fs::Real,freqmin::Float64,freqmax::Float64,
                     cc_step::Int, cc_len::Int, starttime::DateTime,
                     endtime::DateTime; time_norm::Union{Bool,String}=false,
                     to_whiten::Union{Bool,String}=false)

    process_raw!(S,fs)  # demean, detrend, taper, lowpass, downsample
    merge!(S)
    sync!(S,starttime,endtime)
    A, starts, ends = slide(S[1], cc_len, cc_step)
    FFT = process_fft(A, freqmin, freqmax, fs, time_norm=time_norm,
                      to_whiten=to_whiten)
    return F = FFTData(S[1].name, S[1].id, S[1].loc, S[1].fs, S[1].gain, freqmin, freqmax,
                cc_len, cc_step, to_whiten, time_norm, S[1].resp, S[1].misc,
                S[1].notes, hcat(starts,ends), fft)
end

"""

    data_checks!(S::SeisData)

Perform sanity checks on input raw data.
"""
function data_checks!(S::SeisData)
end

"""
    process_raw!(S,fs)

Pre-process raw seismic data.

Checks:
- sample rate is fs
- downsamples data
- checks for gaps in data
- phase-shifts data to begin at 00:00:00.0

# Arguments
- `S::SeisChannel`: SeisData structure.
- `fs::Real`: Sampling rate to downsample `S`.
"""
function process_raw!(S::SeisData, fs::Real)
    demean!(S)        # remove mean from channel
    ungap!(S)         # replace gaps with mean of channel
    detrend!(S)       # remove linear trend from channel
    taper!(S)         # taper channel ends
    lowpass!(S,fs)    # lowpass filter before downsampling
    S = downsample(S,fs) # downsample to lower fs
    phase_shift!(S) # timing offset from sampling period
    return nothing
end
process_raw(S::SeisData, fs::Real) = (U = deepcopy(S);
            process_raw!(U,fs); return U)

"""
    process_fft(A::AbstractArray,freqmin::Real,freqmax::Real,fs::Real;
                time_norm=false,to_whiten=false,corners=corners,
                zerophase=zerophase)

apply 1-bit, filter, whitening
"""
function process_fft(A::AbstractArray,freqmin::Float64,freqmax::Float64,
                     fs::Float64; time_norm::Union{Bool,String}=false,
                     to_whiten::Union{Bool,String}=false,
                     corners::Int=4,
                     zerophase::Bool=true)

    window_samples, N = size(A)

    # pre-process each window
    for ii = 1:N
        ArrayFuncs.demean!(A[:,ii])
        ArrayFuncs.detrend!(A[:,ii])
        taper!(A[:,ii],fs)
        bandpass!(A[:,ii],freqmin,freqmax,fs,corners=corners,
                  zerophase=zerophase)
        ArrayFuncs.demean!(A[:,ii])
    end

    # apply one-bit normalization
    if time_norm == "one_bit"
        A .= sign.(A)
    end

    # take fft along first dimension
    if to_whiten
        FFT = whiten(A,freqmin, freqmax, fs)
    else
        FFT = rfft(A,1)
    end
    return FFT
end


"""
    remove_resp(args)

remove instrument response - will require reading stationXML and extracting poles
and zeros
"""
function remove_resp(args)
end
