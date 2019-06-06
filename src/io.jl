export save_fft, load_fft, load_corr, save_corr

"""
    save_fft(F::FFTData, OUT::String)

Save FFTData `F` to JLD2.
"""
function save_fft(F::FFTData, FFTOUT::String)
    # check if FFT DIR exists
    if isdir(FFTOUT) == false
        mkpath(FFTOUT)
    end

    # create JLD2 file and save FFT
    net,sta,loc,chan = split(F.name,'.')
    filename = joinpath(FFTOUT,"$net.$sta.$chan.jld2")
    file = jldopen(filename, "a+")
    if !(chan in keys(file))
        group = JLD2.Group(file, chan)
        group[F.id] = F
    else
        file[chan][F.id] = F
    end
    close(file)
end

"""

  load_fft(filename,chan)

Loads FFTData for channel `chan` from JLD2 file `filename`.
If day is specified, loads data from day `day`, else
loads data from all days of `chan`.
"""
function load_fft(filename::String,chan::String;day::Union{String,Missing}=missing)
    file = jldopen(filename,"a+")
    if ismissing(day)
        days = keys(file[chan])
        F = FFTData()
        for day in days
            F += file[chan][day]
        end
    else
        F = file[chan][day]
    end
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

"""

  load_corr(filename,chan,day)

Loads CorrData for channel `chan` on day `day` from JLD2 file `filename`.
"""
function load_corr(filename::String,chan::String;
                   day::Union{String,Missing}=missing)
    file = jldopen(filename,"a+")
    if ismissing(day)
        days = keys(file[chan])
        C = CorrData()
        for day in days
            C += file[chan][day]
        end
    else
        C = file[chan][day]
    end
    close(file)
    return C
end
