export process_raw, process_raw!, process_fft, compute_fft
export onebit!, onebit, remove_response!, remove_response, amplitude!, amplitude
export clip, clip!, clamp, clamp!, mute!, mute
import Base: clamp, clamp!
"""
    process_raw!(S,fs)

Pre-process raw seismic data.

- Removes mean from each channel in `S`.
- Detrends each channel in `S`.
- Downsamples data to sampling rate `fs`
- Phase-shifts data to begin at 00:00:00.0

# Arguments
- `S::SeisData`: SeisData structure.
- `fs::Float64`: Sampling rate to downsample `S`.
"""
function process_raw!(S::SeisData, fs::Float64; ϕshift::Bool=true)
    merge!(S)
    ungap!(S)

    for ii = 1:S.n
        SeisIO.detrend!(S[ii])         # remove mean & trend from channel
        if fs ∉ S.fs
            SeisNoise.taper!(S[ii].x,S[ii].fs)         # taper channel ends
            lowpass!(S[ii].x,fs/2,S[ii].fs)    # lowpass filter before downsampling
        end
        resample!(S,chans=ii,fs=fs) # downsample to lower fs
        phase_shift!(S[ii], ϕshift=ϕshift) # timing offset from sampling period
    end
    return nothing
end
process_raw(S::SeisData, fs::Float64;
           ϕshift::Bool=true) = (U = deepcopy(S);
           process_raw!(U,fs, ϕshift=ϕshift); return U)

"""

  compute_fft(R)

Computes windowed rfft of ambient noise data. Returns FFTData structure.

# Arguments
- `R::RawData`: RawData structure
"""
function compute_fft(R::RawData)
    FFT = rfft(R.x,1)
    return FFTData(R.name, R.id,R.loc, R.fs, R.gain, R.freqmin, R.freqmax,
                   R.cc_len, R.cc_step, R.whitened, R.time_norm, R.resp,
                   R.misc, R.notes, R.t, FFT)
end

"""

  mute(A,factor=factor)

Set high amplitudes in array `A` to zero.
Uses median of envelope of `A` to find outliers. 
"""
function mute!(A::AbstractArray;factor::Int=3)
    T = eltype(A)
    envelope = abs.(hilbert(A))
    levels = mean(envelope,dims=1)
    level = factor * median(levels)
    A[envelope .> level] .= T(0)
    return nothing
end
mute(A::AbstractArray;factor::Int=3) = (U = deepcopy(A); mute!(U,factor=factor);
     return U)
mute!(R::RawData;factor::Int=3) = mute!(R.x,factor=factor)
mute(R::RawData;factor::Int=3) = (U = deepcopy(R); mute!(U,factor=factor);
     return U)

"""
  clip(A,factor)

Truncate array A at `factor` times the root mean square of each column.

#Arguments
- `A::AbstractArray`: N-d time series array
- `factor::Real`:
- `f::Function`: Input statistical function (e.g. rms, var, std, mad)
- `dims`: Dimension of `A` to apply clipping (defaults to 1)
"""
function clip!(A::AbstractArray{T,N}, factor::Real; f::Function=rms,dims=1) where {T,N}
    if N == 1
        high = f(A) .* factor
        clamp!(@view(A[:]),-high,high)
    else
        for ii = 1:size(A,2)
            high = f(@view(A[:,ii])) * factor
            clamp!(@view(A[:,ii]),-high,high)
        end
    end
    return nothing
end
clip(A::AbstractArray, factor::Real; f::Function=rms, dims=1) = (U = deepcopy(A);
     clip!(U,factor,f=f,dims=dims);return U)
clip!(R::RawData,factor::Real;f::Function=rms,dims=1) = clip!(R.x,factor,f=f,dims=dims)
clip(R::RawData,factor::Real;f::Function=rms,dims=1) = (U = deepcopy(R);
     clip!(U.x,factor,f=f,dims=dims); return U)
clamp!(R::RawData,val::Real) = clamp!(R.x,-abs(val),abs(val))
clamp!(R::RawData,lo::Real,hi::Real) = clamp!(R.x,lo,hi)
clamp(R::RawData,val::Real) = (U = deepcopy(R); clamp!(U,val);return U)
clamp(R::RawData,lo::Real,hi::Real) = (U = deepcopy(R); clamp!(U,lo,hi);return U)

"""
  amplitude!(R)

Filter raw data based on amplitude.
"""
function amplitude!(R::RawData; max_std::Float64=10.)
    # remove nonzero columns
    zeroind = nonzero(R.x)
    if length(zeroind) == 0
        R = nothing
    elseif size(R.x,2) != length(zeroind)
        R.x = R.x[:,zeroind]
        starts = starts[zeroind]
    end

    # amplitude threshold indices
    stdind = std_threshold(R.x,max_std)
    if length(stdind) == 0
        R = nothing
    elseif size(R.x,2) != length(stdind)
        R.x = R.x[:,stdind]
        starts = starts[stdind]
        ends = ends[stdind]
    end
    return nothing
end
amplitude(R::RawData; max_std::Float64=10.) = (U = deepcopy(R);amplitude!(U,
          max_std=max_std); return U)

"""

  onebit!(R)

One-bit amplitude modification of RawData `R`.
"""
function onebit!(R::RawData)
    R.x .= sign.(R.x)
    return nothing
end
onebit(R::RawData) = (U = deepcopy(R); sign!(U);return U)

"""
  remove_response!(S, stationXML, freqmin, freqmax)

Loads instrument response from stationXML and removes response from `S`.

# Arguments
- `S::SeisData`: SeisData structure.
- `stationXML::String`: Path to stationXML file, e.g. "/path/to/file.xml"
- `freqmin::Float64`: minimum frequency for pre-filtering.
- `freqmax::Float64`: maximum frequency for pre-filtering.
- `np::Int`: number of poles in pre-filter.
- `t_max::Float64`: Length of taper in seconds.
- `wl::Float32`: Water level for instrument response.
"""
function remove_response!(S::SeisData, stationXML::String, freqmin::Float64,
                          freqmax::Float64;np::Int=2, t_max::Float64=20.,
                          wl::Float32=eps(Float32))

    # read response file
    if !isfile(stationXML)
        error("$stationXML does not exist. Instrument response not removed.")
    end

    s,e = string.(start_end(S[1]))
    R = read_sxml(stationXML,s=s,t=e)
    # loop through responses
    Rid = R.id
    Sid = S.id

    # remove trend, taper and prefilter
    SeisIO.detrend!(S)
    SeisIO.taper!(S,t_max=t_max)
    filtfilt!(S,fl=freqmin,fh=freqmax,np=np,rt="Bandpass")

    # remove instrument response for each channel in S
    @inbounds for ii = 1:S.n
        id = S[ii].id
        ind = findfirst(x -> x == id,Rid)
        LOC = R[ind].loc
        GAIN = R[ind].gain
        RESP = R[ind].resp
        UNITS = R[ind].units
        setindex!(S.loc,LOC,ii)
        setindex!(S.gain,GAIN,ii)
        setindex!(S.units,UNITS,ii)
        translate_resp!(S,RESP,chans=ii,wl=wl)
    end
    return nothing
end
remove_response(S::SeisData, stationXML::String, freqmin::Float64,
                          freqmax::Float64;np::Int=2, t_max::Float64=20.,
                          wl::Float32=eps(Float32)) = (U = deepcopy(S);
            remove_response!(U,stationXML,freqmin,freqmax,np=np,t_max=t_max,
                             wl=wl); return U)

"""
nonzero(A)

Find indices of all nonzero columns in array `A`.
"""
function nonzero(A::AbstractArray)
    Nrows, Ncols = size(A)
    ind = Int64[]
    sizehint!(ind,Ncols)
    for ii = 1:Ncols
        for jj = 1:Nrows
            if !iszero(A[jj,ii])
                append!(ind,ii)
                break
            end
        end
    end
    return ind
end

"""
    compute_fft(S, freqmin, freqmax, fs, cc_step, cc_len;
                time_norm=false, to_whiten=false, max_std=5.)
Computes windowed fft of ambient noise data.
# Arguments
- `S::SeisData`: SeisData structure.
- `freqmin::Float64`: minimum frequency for filtering/whitening.
- `freqmax::Float64`: maximum frequency for filtering/whitening.
- `fs::Float64`: Sampling rate to downsample `S`.
- `cc_step::Int`: time, in seconds, between successive cross-correlation windows.
- `cc_len::Int`: length of noise data window, in seconds, to cross-correlate.
- `time_norm::Union{Bool,String}`: time domain normalization to perform.
- `to_whiten::Bool`: Apply whitening in frequency domain.
- `max_std::Float64=5.`: Number of standard deviations above mean to reject windowed data.
"""
function compute_fft(S::SeisData,freqmin::Float64,freqmax::Float64,fs::Float64,
                     cc_step::Int, cc_len::Int;
                     time_norm::Union{Bool,String}=false,
                     to_whiten::Bool=false,
                     max_std::Float64=5.,
                     ϕshift::Bool=true)

    # sync!(S,s=starttime,t=endtime)
    merge!(S)
    ungap!(S)
    process_raw!(S,fs,ϕshift=ϕshift)  # demean, detrend, taper, lowpass, downsample

    # subset by time
    starttime, endtime = u2d.(nearest_start_end(S[1],cc_len, cc_step))
    # check if waveform length is < cc_len
    if Int(floor((endtime - starttime).value / 1000)) < cc_len
        return nothing
    end

    sync!(S,s=starttime,t=endtime)
    A, starts, ends = slide(S[1], cc_len, cc_step)

    # remove nonzero columns
    zeroind = nonzero(A)
    if length(zeroind) == 0
        return nothing
    elseif size(A,2) != length(zeroind)
        A = A[:,zeroind]
        starts = starts[zeroind]
        ends = ends[zeroind]
    end

    # amplitude threshold indices
    stdind = std_threshold(A,max_std)
    if length(stdind) == 0
        return nothing
    elseif size(A,2) != length(stdind)
        A = A[:,stdind]
        starts = starts[stdind]
        ends = ends[stdind]
    end

    FFT = process_fft(A, freqmin, freqmax, fs, time_norm=time_norm,
                      to_whiten=to_whiten)
    return F = FFTData(S[1].id, Dates.format(u2d(starts[1]),"Y-mm-dd"),
                       S[1].loc, S[1].fs, S[1].gain, freqmin, freqmax,
                       cc_len, cc_step, to_whiten, time_norm, S[1].resp,
                       S[1].misc, S[1].notes, starts, FFT)
end
