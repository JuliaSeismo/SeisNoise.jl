export process_raw, process_raw!, process_fft, compute_fft, save_fft

"""
    compute_fft()


Computes windowed fft of ambient noise data.

Cross-correlates data using either cross-correlation, deconvolution or
cross-coherence. Saves cross-correlations in JLD2 data set.

TO DO:
    - load in data [x]
    - process_raw [x]
    - check start/end times [x]
    - chop into matrix [x]
    - normalize time / freq domain [x]
    - take fft [x]
    - get parameters for each window (amplitude, mad) []
    - save fft and parameters to JLD2 [x]


:type maxlag: int
:param maxlag: maximum lag, in seconds, in cross-correlation
:type fs: Float64
:param fs: Frequency to which waveforms in stream are downsampled
:type freqmin: float
:param freqmin: minimun frequency for whitening
:type freqmax: float
:param freqmax: maximum frequency for whitening
:type cc_step: Int
:param cc_step: time, in seconds, between success cross-correlation windows
:type cc_len: Int
:param cc_len: length of noise data window, in seconds, to cross-correlate
"""
function compute_fft(S::SeisData,fs::Float64,freqmin::Float64,freqmax::Float64,
                     cc_step::Int, cc_len::Int;
                     time_norm::Union{Bool,String}=false,
                     to_whiten::Union{Bool,String}=false)

    # sync!(S,s=starttime,t=endtime)
    process_raw!(S,fs)  # demean, detrend, taper, lowpass, downsample
    merge!(S)
    starttime, endtime = u2d.(nearest_start_end(S[1],cc_len, cc_step))
    sync!(S,s=starttime,t=endtime)
    A, starts, ends = slide(S[1], cc_len, cc_step)
    FFT = process_fft(A, freqmin, freqmax, fs, time_norm=time_norm,
                      to_whiten=to_whiten)
    return F = FFTData(S[1].name, Dates.format(u2d(starts[1]),"Y-mm-dd"),
                       S[1].loc, S[1].fs, S[1].gain, freqmin, freqmax,
                       cc_len, cc_step, to_whiten, time_norm, S[1].resp,
                       S[1].misc, S[1].notes, starts, FFT)
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
- `fs::Float64`: Sampling rate to downsample `S`.
"""
function process_raw!(S::SeisData, fs::Float64)
    demean!(S)        # remove mean from channel
    ungap!(S)         # replace gaps with mean of channel
    detrend!(S)       # remove linear trend from channel
    taper!(S)         # taper channel ends
    lowpass!(S,fs/2)    # lowpass filter before downsampling
    S = downsample(S,fs) # downsample to lower fs
    phase_shift!(S) # timing offset from sampling period
    return nothing
end
process_raw(S::SeisData, fs::Float64) = (U = deepcopy(S);
            process_raw!(U,fs); return U)

"""
    process_fft(A::AbstractArray,freqmin::Float64,freqmax::Float64,fs::Float64;
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

"""
    save_fft(F::FFTData, OUT::String)

Save FFTData `F` to JLD2.
"""
function save_fft(F::FFTData, FFTOUT::String)
    # check if FFT DIR exists
    if isdir(FFTOUT) == false
        mkpath(FFTOUT)
    end

    # create JLD2 file and save mseed
    net,sta,loc,chan = split(F.name,'.')
    filename = joinpath(FFTOUT,"$net.$sta.jld2")
    file = jldopen(filename, "a+")
    if !(chan in keys(file))
        group = JLD2.Group(file, chan)
        group[F.id] = F
    else
        file[chan][F.id] = F
    end
    close(file)
end
