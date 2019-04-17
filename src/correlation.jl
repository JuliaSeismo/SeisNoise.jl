# cross-correlation module
export clean_up!, clean_up, correlate, compute_cc, next_fast_len
import ArrayFuncs
"""
    clean_up!(A,fs,freqmin,freqmax)

Demean, detrend, taper and filter time series.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
"""
function clean_up!(A::AbstractArray, fs::Real, freqmin::Real, freqmax::Real;
                   corners::Int=4, zerophase::Bool=true)
    ArrayFuncs.demean!(A)
    ArrayFuncs.detrend!(A)
    ArrayFuncs.taper!(A,fs)
    ArrayFuncs.bandpass!(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
    return nothing
end
clean_up(A::AbstractArray, fs::Real, freqmin::Real, freqmax::Real;
         corners::Int=4, zerophase::Bool=false) =
                 (U = deepcopy(A); clean_up!(U,fs,freqmin,freqmax,
                  corners=corners, zerophase=zerophase); return U)

"""
    correlate(fft1, fft2, N, maxlag, corr_type='cross-correlate')

Cross-correlation of two ffts.


"""
function correlate(fft1::AbstractArray, fft2::AbstractArray, N::Int,
                   maxlag::Int;
                   smoothing_half_win::Int=20,
                   corr_type::String="cross-correlation")

    corrF = fft1 .* conj(fft2)
    if corr_type == "deconv"
        corrF ./= (smooth(abs.(fft2).^2, half_win=smoothing_half_win) .+
                   0.01 .* mean(smooth(abs.(fft2).^2, half_win=smoothing_half_win),dims=2))
    elseif corr_type == "coherence"
        corrF ./= smooth(abs.(fft1),half_win=smoothing_half_win) .+
                   0.01 .* mean(smooth(abs.(fft1), half_win=smoothing_half_win),dims=2)
        corrF ./= smooth(abs.(fft2),half_win=smoothing_half_win) .+
                   0.01 .* mean(smooth(abs.(fft2), half_win=smoothing_half_win),dims=2)
    end

    # take inverse fft
    corrT = irfft(corrF,N)
    corrT = fftshift(corrT)

    # return corr[-maxlag:maxlag]
    t = range(-Int(N/2) + 1, Int(N/2) - 1)
    ind = findall(x -> abs(x) <= maxlag,t)
    corrT = corrT[ind,:]
end

"""
    compute_cc(FFT1::FFTData, FFT2::FFTData, N::Int, maxlag::Float64,
               comp::String; smoothing_half_win::Int=20,
               corr_type::String="cross-correlation" )

"""
function compute_cc(FFT1::FFTData, FFT2::FFTData, N::Int, maxlag::Float64,
                    comp::String;
                    smoothing_half_win::Int=20,
                    corr_type::String="cross-correlation")
    # get intersect of dates; return nothing if no intersect
    inter = intersect(FFT1.t,FFT2.t)
    if length(inter) == 0
        return nothing
    end

    ind1 = findall(x -> x ∈ inter, FFT1.t)
    ind2 = findall(x -> x ∈ inter, FFT2.t)

    corr = correlate(FFT1.fft[:,ind1], FFT2.fft[:,ind2], N,
                     convert(Int,round(maxlag * FFT1.fs)),
                     corr_type=corr_type)
    rotated = false

    return CorrData(FFT1, FFT2, comp, rotated, corr_type,
                    maxlag, inter, corr)

end

"""
    next_fast_len(N::Real)

Return next fast length for fft with FFTW.
"""
function next_fast_len(N::Real)
    return nextprod([2,3,5],N)
end

"""
    save_corr(C::CorrData, OUT::String)

Save CorrData `C` to JLD2.
"""
function save_corr(C::CorrData, CORROUT::String)
    # check if FFT DIR exists
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end

    # create JLD2 file and save correlation
    net1,sta1,loc1,chan1,net2,sta2,loc1,chan2 = split(C.name,'.')
    filename = joinpath(CORROUT,"$net1.$sta1.$net2.$sta2.jld2")
    file = jldopen(filename, "a+")
    if !(chan in keys(file))
        group = JLD2.Group(file, chan)
        group[C.id] = C
    else
        file[chan][C.id] = C
    end
    close(file)
end
