export raw_data_available

function raw_data_available(ROOT::String; format::String="mseed")
    files = glob("*"*format,ROOT)
    N = length(files)
    NET = Array{String,1}(undef,N)
    STA = Array{String,1}(undef,N)
    LOC = Array{String,1}(undef,N)
    CHAN = Array{String,1}(undef,N)
    STARTTIME = Array{DateTime,1}(undef,N)
    ENDTIME = Array{DateTime,1}(undef,N)

    for ii = 1:N
        path, mseed = splitdir(files[ii])
        chan, LOC[ii], starttime, endtime, _ = string.(split(mseed,'.'))
        path, CHAN[ii] = splitdir(path)
        path, STA[ii] = splitdir(path)
        path, NET[ii] = splitdir(path)
        STARTTIME[ii] = DateTime(starttime,"yyyymmddTHHMMSSZ")
        ENDTIME[ii] = DateTime(endtime,"yyyymmddTHHMMSSZ")
    end

    df = DataFrame(NET = NET, STA = STA, LOC = LOC, CHAN = CHAN,
               STARTTIME = STARTTIME, ENDTIME = ENDTIME, FILE = files)
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
