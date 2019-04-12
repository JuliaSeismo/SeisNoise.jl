export get_files, get_info, raw_data_available

function get_files(ROOT::String; format::String="mseed")
    # add files to path
    paths = []
    for (root, dirs, files) in walkdir(ROOT)
        for file in files
            file = joinpath(root,file)
            push!(paths,file)
        end
    end

    # check file format
    filter!(x-> occursin(Regex(".$format"),x), paths)
    return paths
end

function get_info(file::String)
    path, mseed = splitdir(file)
    chan, loc, starttime, endtime, _ = split(mseed,'.')
    path, chan = splitdir(path)
    path, sta = splitdir(path)
    path, net = splitdir(path)
    starttime = DateTime(starttime,"yyyymmddTHHMMSSZ")
    endtime = DateTime(endtime,"yyyymmddTHHMMSSZ")
    return Dict(:NET => net,:STA => sta, :LOC => loc, :CHAN => chan,
                :STARTTIME => starttime, :ENDTIME => endtime,
                :FILE => file)
end

function raw_data_available(ROOT::String; format::String="mseed")
    files = get_files(ROOT)
    df = DataFrame(NET = String[], STA = String[], LOC = String[], CHAN = String[],
               STARTTIME = DateTime[], ENDTIME = DateTime[], FILE = String[])

    # daily files to data frame
    for ii = 1:length(files)
        push!(df,get_info(files[ii]))
    end
    return df
end

function mseed2jld2(df::DataFrame, DIR::String)

    # check if base dir is in file
    if isdir(DIR) == false
        mkpath(DIR)
    end

    # group stations by network and station
    for station in groupby(df, [:NET, :STA, :CHAN])
        net = unique(station[:NET])[1]
        sta = unique(station[:STA])[1]
        chan = unique(station[:CHAN])[1]
        NET = joinpath(DIR,net)

        if isdir(NET) == false
            mkpath(NET)
        end

        # create JLD2 file and save mseed
        filename = joinpath(NET, "$net.$sta.jld2")
        file = jldopen(filename, "a+")
        group = JLD2.Group(file, chan)
        println("Adding station $net.$sta.$chan $(Dates.now())")

        for day in eachrow(station)
            st = readmseed(day[:FILE])
            group[Dates.format(day[:STARTTIME],"yyyy.mm.dd")] = st
        end

        close(file)
    end
end

function readjld2(filename::String, chan::String, starttime::DateTime,
                  endtime::DateTime)
    file = jldopen(filename,"r")
    days = keys(file[chan])
    dates = [DateTime(z,"yyyy.mm.dd") for z in days]
    ind = findall(x -> (x >= starttime) & (x <= endtime), dates)
    days = days[ind]
    data = Array{SeisData}(undef,length(days))
    for (ii,day) in enumerate(days)
        data[ii] = file[chan][day]
    end
    close(file)
    return data
end

function date_range(dates::AbstractArray, starttime::DateTime, endtime::DateTime)
    return [d for d in dates if (d >= starttime) & (d <= endtime)]
end
