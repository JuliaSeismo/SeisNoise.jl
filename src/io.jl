export save_fft, load_fft, load_corr, save_corr

"""
    save_fft(F::FFTData, OUT::String)

Save FFTData `F` to JLD2.
"""
function save_fft(F::FFTData, FFTOUT::String)
    # check if FFT DIR exists
    FFTOUT = expanduser(FFTOUT)
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

  load_fft(filename,chan,day=day)

Loads FFTData for channel `chan` from JLD2 file `filename`.
If day is specified, loads data from day `day`, else
loads data from all days of `chan`.
"""
function load_fft(filename::String,chan::String;day::Union{String,Missing}=missing)
    file = jldopen(filename,"r")
    if ismissing(day)
        days = keys(file[chan])
        num_days  = length(days)

        # get sizes of all corrs
        sizes = Array{Int}(undef,num_days,2)

        for ii = 1:num_days
            sizes[ii,:]=collect(size(file[chan][days[ii]].fft))
        end

        # create empty arrays and get indicies
        M = sizes[1]
        N = sum(sizes[:,2])
        ffts = Array{eltype(file[chan][days[1]].fft)}(undef,M,N)
        t = Array{eltype(file[chan][days[1]].t)}(undef,N)
        starts = cumsum(vcat(1 ,sizes[:,2]))[1:end-1]
        ends = cumsum(sizes[:,2])
        for ii = 1:num_days
            ffts[:,starts[ii]:ends[ii]] = file[chan][days[ii]].fft
            t[starts[ii]:ends[ii]] = file[chan][days[ii]].t
        end

        # add data to FFTData
        F = file[chan][days[1]]
        F.t = t
        F.fft = ffts
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
    CORROUT = expanduser(CORROUT)
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

  load_corr(filename,chan,day=day)

Loads CorrData for channel `chan` on day `day` from JLD2 file `filename`.
"""
function load_corr(filename::String,chan::String;
                   day::Union{String,Missing}=missing)
    file = jldopen(filename,"r")
    if ismissing(day) # load all files
        days = keys(file[chan])
        num_days  = length(days)

        # get sizes of all corrs
        sizes = Array{Int}(undef,num_days,2)

        for ii = 1:num_days
            sizes[ii,:]=collect(size(file[chan][days[ii]].corr))
        end

        # create empty arrays and get indicies
        M = sizes[1]
        N = sum(sizes[:,2])
        corrs = Array{eltype(file[chan][days[1]].corr)}(undef,M,N)
        t = Array{eltype(file[chan][days[1]].t)}(undef,N)
        starts = cumsum(vcat(1 ,sizes[:,2]))[1:end-1]
        ends = cumsum(sizes[:,2])
        for ii = 1:num_days
            corrs[:,starts[ii]:ends[ii]] = file[chan][days[ii]].corr
            t[starts[ii]:ends[ii]] = file[chan][days[ii]].t
        end

        # add data to CorrData
        C = file[chan][days[1]]
        C.t = t
        C.corr = corrs
    else # return single day
        C = file[chan][day]
    end
    close(file)
    return C
end
