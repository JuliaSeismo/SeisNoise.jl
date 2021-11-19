export process_raw, process_raw!, process_fft, rfft
export onebit!, onebit, remove_response!, remove_response, remove_amp!, remove_amp
export clip, clip!, clamp, clamp!, mute!, mute, phase
import Base: clamp, clamp!
import FFTW: rfft
"""
    process_raw!(S,fs)

Pre-process raw seismic data.

- Removes mean from each channel in `S`.
- Detrends each channel in `S`.
- Downsamples data to sampling rate `fs`
- Phase-shifts data to begin at 00:00:00.0

# Arguments
- `S::SeisData`: SeisData structure.
- `fs::Real`: Sampling rate to downsample `S`.
"""
function process_raw!(S::SeisData, fs::Real; ϕshift::Bool=true)
    merge!(S)
    ungap!(S)
	detrend!(S)         # remove mean & trend from channel
	taper!(S)                      # taper channel ends
	if fs ∉ S.fs
		filtfilt!(S,fh=Float64(fs/2),rt="Lowpass")    # lowpass filter before downsampling
	end
	resample!(S,fs=Float64(fs)) # downsample to lower fs
	taper!(S)
    phase_shift!(S, ϕshift=ϕshift) # timing offset from sampling period
    return nothing
end
process_raw(S::SeisData, fs::Real;
           ϕshift::Bool=true) = (U = deepcopy(S);
           process_raw!(U,fs, ϕshift=ϕshift); return U)

"""
   process_raw!(C,fs)

Pre-process raw SeisChannel.

- Removes mean from data in `C`.
- Detrends data in `C`.
- Downsamples data to sampling rate `fs`
- Phase-shifts data to begin at 00:00:00.0

# Arguments
- `C::SeisChannel`: SeisChannel structure.
- `fs::Real`: Sampling rate to downsample `C`.
"""
function process_raw!(C::SeisChannel, fs::Real; ϕshift::Bool=true)
	ungap!(C)
    detrend!(C)         # remove mean & trend from channel
    taper!(C)         # taper channel ends
    if fs != C.fs
    	filtfilt!(C,fh=Float64(fs/2),rt="Lowpass")    # lowpass filter before downsampling
	end
	resample!(C,Float64(fs)) # downsample to lower fs
	taper!(C)
	phase_shift!(C, ϕshift=ϕshift) # timing offset from sampling period
	return nothing
end
process_raw(C::SeisChannel, fs::Real;
        ϕshift::Bool=true) = (U = deepcopy(C);
        process_raw!(U,fs, ϕshift=ϕshift); return U)

"""

 rfft(R)

Computes windowed rfft of ambient noise data. Returns FFTData structure.

Overloads the rfft function from FFTW.

# Arguments
- `R::RawData`: RawData structure
"""
function rfft(R::RawData,dims::Int=1)
    FFT = rfft(R.x,dims)
    return FFTData(R.name, R.id,R.loc, R.fs, R.gain, R.freqmin, R.freqmax,
                 R.cc_len, R.cc_step, R.whitened, R.time_norm, R.resp,
                 R.misc, R.notes, R.t, FFT)
end

"""
    phase(A::AbstractArray)

Extract instantaneous phase from signal A.

For time series `A`, its analytic representation ``S = A + H(A)``, where
``H(A)`` is the Hilbert transform of `A`. The instantaneous phase ``e^{iθ}``
of `A` is given by dividing ``S`` by its modulus: ``e^{iθ} = \\frac{S}{|S|}``
For more information on Phase Cross-Correlation, see:
[Ventosa et al., 2019](https://pubs.geoscienceworld.org/ssa/srl/article-standard/570273/towards-the-processing-of-large-data-volumes-with).
"""
function phase(A::AbstractArray)
    # the analytic signal 
    s =  analytic(A)
    return s ./ abs.(s)
end

function analytic(A::AbstractArray)
    # the analytic signal 
    T = real(eltype(A))
    return hilberttransform(A) .* Complex(T(0),T(1)) .+ A
end

function hilberttransform(A::AbstractArray)
    Nrows = size(A,1)
    T = real(eltype(A))
    f = fft(A,1)
    f[1,:] .*= Complex(T(0),T(0))
    if iseven(Nrows)
        f[2:Nrows÷2 + Nrows % 2,:] .*= Complex(T(0),T(-1))
        f[Nrows÷2 + Nrows % 2 + 1,:] .*= Complex(T(0),T(0))
        f[Nrows÷2 + Nrows % 2 + 2: end,:] .*= Complex(T(0),T(1))
    else
        f[2:Nrows÷2 + Nrows % 2,:] .*= Complex(T(0),T(-1))
        f[Nrows÷2 + Nrows % 2 + 1 : end,:] .*= Complex(T(0),T(1))
    end
    return ifft(f,1)
end

"""

 phase(R)

Computes windowed analytic signal of ambient noise data. Returns FFTData structure.


# Arguments
- `R::RawData`: RawData structure
"""
function phase(R::RawData)
    FFT = phase(R.x)
    return FFTData(R.name, R.id,R.loc, R.fs, R.gain, R.freqmin, R.freqmax,
                 R.cc_len, R.cc_step, R.whitened, R.time_norm, R.resp,
                 R.misc, R.notes, R.t, FFT)
end

isweird(x) =  isnan(x) .| isinf(x)

"""

  mute(A,factor)

Set high amplitudes in array `A` to zero.
Uses median of envelope of `A` to find outliers.
"""
function mute!(A::AbstractArray,factor::Real=3)
    T = eltype(A)
    envelope = abs.(hilbert(A))
    levels = mean(envelope,dims=1)
    level = factor .* median(levels)
    A[envelope .> level] .= T(0)
    return nothing
end
mute(A::AbstractArray,factor::Real=3) = (U = deepcopy(A); mute!(U,factor);
     return U)
mute!(R::RawData,factor::Real=3) = mute!(R.x,factor)
mute(R::RawData,factor::Real=3) = (U = deepcopy(R); mute!(U,factor);
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
function clip!(A::AbstractArray{T,N}, factor::Real; f::Function=std,dims=1) where {T,N}
    if N == 1
        high = f(A) .* factor
        clamp!(@view(A[:]),-high,high)
    else
        high = f(@view(A[:,:]),dims=dims) .* factor
        for ii = 1:size(A,2)
            clamp!(@view(A[:,ii]),-high[ii],high[ii])
        end
    end
    return nothing
end
clip(A::AbstractArray, factor::Real; f::Function=std, dims=1) = (U = deepcopy(A);
     clip!(U,factor,f=f,dims=dims);return U)
clip!(R::RawData,factor::Real;f::Function=std,dims=1) = clip!(R.x,factor,f=f,dims=dims)
clip(R::RawData,factor::Real;f::Function=std,dims=1) = (U = deepcopy(R);
     clip!(U.x,factor,f=f,dims=dims); return U)
clamp!(R::RawData,val::Real) = clamp!(R.x,-abs(val),abs(val))
clamp!(R::RawData,lo::Real,hi::Real) = lo < hi ? clamp!(R.x,lo,hi) : throw(ArgumentError("lo value $lo must be less than hi value $hi."))
clamp(R::RawData,val::Real) = (U = deepcopy(R); clamp!(U,val);return U)
clamp(R::RawData,lo::Real,hi::Real) = (U = deepcopy(R); clamp!(U,lo,hi);return U)

"""
  remove_amp!(R)

Filter raw data based on amplitude.
"""
function remove_amp!(R::RawData; max_std::Real=10.)
    # remove nonzero columns
    zeroind = nonzero(R.x)
    if length(zeroind) == 0
        throw(ErrorException("All values in Rawdata == 0"))
    elseif size(R.x,2) != length(zeroind)
        R.x = R.x[:,zeroind]
        R.t = R.t[zeroind]
    end

    # amplitude threshold indices
    stdind = std_threshold(R.x,max_std)
    if length(stdind) == 0
        throw(ErrorException("All columns in Rawdata contain values > max_std"))
    elseif size(R.x,2) != length(stdind)
        R.x = R.x[:,stdind]
        R.t = R.t[stdind]
    end
    return nothing
end
remove_amp(R::RawData; max_std::Real=10.) = (U = deepcopy(R);remove_amp!(U,
          max_std=max_std); return U)

"""

  onebit!(R)

One-bit amplitude modification of RawData `R`.
"""
function onebit!(R::RawData)
    R.x .= sign.(R.x)
    return nothing
end
onebit(R::RawData) = (U = deepcopy(R); onebit!(U);return U)

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
