# cross-correlation module
export clean_up!, clean_up, correlate, compute_cc, next_fast_len, save_corr, load_fft
export correlate_parallel, generate_pairs
import ArrayFuncs
import Statistics: mean
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
    t = range(-Int(N/2) + 1, stop=Int(N/2) - 1)
    ind = findall(x -> abs(x) <= maxlag,t)
    corrT = corrT[ind,:]
end

"""
    compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;
               smoothing_half_win::Int=20,
               corr_type::String="cross-correlation" )

"""
function compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;
                    smoothing_half_win::Int=20,
                    corr_type::String="cross-correlation")

    N = convert(Int,round(FFT1.cc_len * FFT1.fs)) # number of data points
    comp = FFT1.name[end] * FFT2.name[end]
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

  load_fft(filename,chan,day)

Loads FFTData for channel `chan` on day `day` from JLD2 file `filename`.
"""
function load_fft(filename::String,chan::String,day::String)
    file = jldopen(filename,"a+")
    F = file[chan][day]
    close(file)
    return F
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
    filename = joinpath(CORROUT,"$net1.$sta1.$chan1.$net2.$sta2.$chan2.jld2")
    file = jldopen(filename, "a+")
    if !(C.comp in keys(file))
        group = JLD2.Group(file, C.comp)
        group[C.id] = C
    else
        file[C.comp][C.id] = C
    end
    close(file)
end

function correlate_parallel(pair,maxlag,smoothing_half_win,corr_type,FFTOUT,CORROUT)
    nsc1, nsc2 = pair[1],pair[2]
    net1,sta1,chan1 = string.(split(nsc1,"."))
    net2,sta2,chan2 = string.(split(nsc2,"."))
    filepath1 = joinpath(FFTOUT,nsc1*".jld2")
    filepath2 = joinpath(FFTOUT,nsc2*".jld2")
    println("Correlating $nsc1 with $nsc2")

    # autocorrelation vs cross-correlation
    if filepath1 == filepath2
        file1 = jldopen(filepath1,"a+")
        days = keys(file1[chan1])
        close(file1)
        for day in days
            FFT1 = load_fft(filepath1,chan1,day)
            C = compute_cc(FFT1,FFT1,maxlag,
                           smoothing_half_win=smoothing_half_win,
                           corr_type=corr_type)
            save_corr(C,CORROUT)
        end
    else # cross-correlations
        file1 = jldopen(filepath1,"a+")
        file2 = jldopen(filepath2,"a+")
        days1 = keys(file1[chan1])
        days2 = keys(file2[chan2])
        close(file1)
        close(file2)
        days = intersect(days1,days2)
        for day in days
            FFT1 = load_fft(filepath1,chan1,day)
            FFT2 = load_fft(filepath2,chan2,day)
            C = compute_cc(FFT1,FFT2,maxlag,
                           smoothing_half_win=smoothing_half_win,
                           corr_type=corr_type)
            save_corr(C,CORROUT)
        end
    end
    return nothing
end

function generate_pairs(files::AbstractArray)
    N = length(files)
    num_pairs = convert(Int,round(N * (N-1) / 2 + N))
    pairs = Array{Array{String,1},1}(undef,num_pairs)
    count = 0
    for ii = 1:length(files)
        for jj = ii:length(files)
            count += 1
            pairs[count] = [files[ii], files[jj]]
        end
    end
    return pairs
end
