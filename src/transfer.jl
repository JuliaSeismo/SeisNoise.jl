export scedctransfer

function df_subset(df::DataFrame,col::String,colsymbol::Symbol)
        col = regex_helper(col)
        ind = findall(occursin.(col,df[!,colsymbol]))
        df = df[ind,:]
	return df
end

function df_subset(df::DataFrame,col::Nothing,colsymbol::Symbol)
	 return df
end

function df_remove(df::DataFrame,col::String,colsymbol::Symbol)
        ind = findall(.!occursin.(col,df[!,colsymbol]))
        df = df[ind,:]
	return df
end

function regex_helper(reg::String)
    if reg == '*'
        # pass for all
    elseif occursin('?',reg)
        reg = replace(reg, '?' => '.')
        reg = Regex(reg)
    elseif occursin('*',reg)
        if reg[end] == '*'
                reg = '^' * strip(reg,'*')
        elseif reg[1] == '*'
                reg = strip(reg,'*') * '$'
        end
            reg = Regex(reg)
    end
    return reg
end

function s3_file_map(filein::String,fileout::String,aws::Dict)
    s3_get_file(aws, "scedc-pds", filein, fileout)
    println("Downloading file: $filein       \r")
end

function filenamecleaner(filename::String)
	filename = strip(filename,'b')
	return replace(filename,"'"=>"")
end

"""

  yearday2Date(year,day)

Convert year and day of year date format to `Date` object.
"""
function yearday2Date(year::String,day::String)
    day = parse(Int,day)
    return Date(year) + Day(day -1)
end

"""

  Date2yearday(d)

Convert `Date` object to yearday string, e.g. 2017354.
"""
function Date2yearday(d::Date)
    days = (d - Date(Year(d))).value + 1
    n = ndigits(days)
    return string(Year(d).value) * ('0' ^ (3 - n)) * string(days)
end

function indexpath(d::Date)
    days = (d - Date(Year(d))).value + 1
    n = ndigits(days)
    outstring = "continuous_waveforms/index/csv/year="
    outstring *= string(Year(d).value) * "/year_doy="
    outstring *= string(Year(d).value) * '_' * ('0' ^ (3 - n)) * string(days)
    outstring *= "/index.csv"
    return outstring
end

function scedcpath(d::Date)
    days = (d - Date(Year(d))).value + 1
    n = ndigits(days)
    outstring = "continuous_waveforms/" * string(Year(d).value) *'/'
    outstring *= string(Year(d).value) * '_' * ('0' ^ (3 - n)) * string(days)
    return outstring
end

"""

  scedctransfer(OUTDIR,startdate)

Transfer data from a single day in the SCEDC open data set from S3 to EC2.

# Arguments
- `OUTDIR::String`: The output directory.
- `startdate::Date`: The start day of the download.
- `enddate::Date`: The (optional) end day of the download.
- `network::String`: Network to download from. If network = "*" or is unspecified,
                       data is downloaded from all available networks.
- `station::String`: Station to download, e.g. "RFO". If station = "*" or is unspecified,
                       data is downloaded from all available stations.
- `channel::String`: Channels to download, e.g. "HH*". If channel = "*" or is unspecified,
                       data is downloaded from all available channels.
- `location::String`: Locations to download, e.g. "00". If channel = "*" or is unspecified,
                       data is downloaded from all available locations. NOTE: most files do
                       not have a location.
- `minlatitude::Float64`: Minimum latitude in data search.
- `maxlatitude::Float64`: Maximum latitude in data search.
- `minlongitude::Float64`: Minimum longitude in data search.
- `maxlongitude::Float64`: Maximum longitude in data search.
"""
function scedctransfer(OUTDIR::String,
                  startdate::Date;
			  	  enddate::Union{Date,Nothing}=nothing,
                  network::Union{String,Nothing}=nothing,
                  station::Union{String,Nothing}=nothing,
                  location::Union{String,Nothing}=nothing,
                  channel::Union{String,Nothing}=nothing,
                  minlatitude::Union{Float64,Nothing}=nothing,
                  maxlatitude::Union{Float64,Nothing}=nothing,
                  minlongitude::Union{Float64,Nothing}=nothing,
                  maxlongitude::Union{Float64,Nothing}=nothing,
                  latitude::Union{Float64,Nothing}=nothing,
                  longitude::Union{Float64,Nothing}=nothing,
                  stationXML::Bool=false)

    tstart = now()
    # connect to S3
    println("Connecting to AWS...      $(now())")
    aws = aws_config(region = "us-west-2")
	@eval @everywhere aws=$aws
    firstdate = Date(2000,1,1)

	if isnothing(enddate)
		enddate = startdate
	end

    # check dates
    if startdate > now()
        @warn("Date must be earlier than today. Aborting download.")
        return nothing
    end

    if enddate < firstdate
        @warn("End date must be later than $firstdate. Aborting download.")
		return nothing
    end

    # get filedf
	date_range = startdate:Day(1):enddate
	paths = indexpath.(date_range)

	println("Starting Download...      $(now())")
	println("Using $(nworkers()) cores...")

	for ii in eachindex(paths)
	    filedf = CSV.read(IOBuffer(s3_get(aws,"scedc-pds",paths[ii])))

	    # subset dataframe
	    # filter by lat/lon
	    if !isnothing(minlatitude)
	        filedf = filedf[filedf[:lat] .> minlatitude,:]
	    end

	    if !isnothing(maxlatitude)
	        filedf = filedf[filedf[:lat] .< maxlatitude,:]
	    end

	    if !isnothing(minlongitude)
	        filedf = filedf[filedf[:lon] .> minlongitude,:]
	    end

	    if !isnothing(maxlongitude)
	        filedf = filedf[filedf[:lon] .< maxlongitude,:]
	    end

	    # filter stations
	    filedf = df_subset(filedf,network,:net)
	    filedf = df_subset(filedf,station,:sta)
	    filedf = df_subset(filedf,channel,:seedchan)
	    filedf = df_subset(filedf,location,:loc)

	    # return if nothing in dataframe
	    if size(filedf,1) == 0
	        @warn("No data available for request $(date_range[ii])!")
			continue
	    end

	    # create directory for instrument responses
	    if stationXML
	        println("Downloading stationXML... $(now())")
	        XMLDIR = joinpath(OUTDIR,"FDSNstationXML")
	        mkpath(XMLDIR)
	        networks = filedf[!,:net]
	        stations = filedf[!,:sta]
	        xmlfiles = networks .* "_" .* stations .* ".xml"
	        ind = indexin(unique(xmlfiles), xmlfiles)
	        networks, stations, xmlfiles = networks[ind], stations[ind], xmlfiles[ind]
	        xml_in = [joinpath("FDSNstationXML",networks[ii],xmlfiles[ii]) for ii = 1:length(xmlfiles)]
	        xml_out = [joinpath(XMLDIR,networks[ii],xmlfiles[ii]) for ii = 1:length(xmlfiles)]
	        xml_dir = unique([dirname(f) for f in xml_out])
	        for ii = 1:length(xml_dir)
	            if !isdir(xml_dir[ii])
	                mkpath(xml_dir[ii])
	            end
	        end


	        # check if requested channels have an instrument response
	        stations2remove = []
	        for ii = 1:length(xml_in)
	            if s3_exists(aws,"scedc-pds",xml_in[ii])
	                s3_get_file(aws,"scedc-pds",xml_in[ii],xml_out[ii])
	            else
	                push!(stations2remove,xmlfiles[ii])
	            end
	        end

	        stations2remove = [split(s[1:end-4],'_')[2] for s in stations2remove]
	        for ii = 1:length(stations2remove)
	            df = df_remove!(filedf,stations2remove[ii],:sta)
	        end
	    end

	    # return if nothing in dataframe
	    if size(filedf,1) == 0
			@warn("No data available for request $(date_range[ii])!")
			continue
	    end

	    # query files
	    dpath = scedcpath(date_range[ii])
	    filenames = filenamecleaner.(filedf[!,:ms_filename])
	    files2download = [joinpath(dpath,f) for f in filenames]

	    # create directory structure
	    OUTDIR = expanduser(OUTDIR)
	    out_files = [joinpath(OUTDIR,dpath,f) for f in filenames]
	    file_dir = unique([dirname(f) for f in out_files])
	    for ii = 1:length(file_dir)
	        if !isdir(file_dir[ii])
	            mkpath(file_dir[ii])
	        end
	    end

	    # download files
	    @eval @everywhere files2download=$files2download
	    @eval @everywhere out_files=$out_files
	    pmap(s3_file_map,files2download,out_files,fill(aws,length(out_files)))
	end
    println("Download Complete!        $(now())          ")
    tend = now()
    println("Download took $(Dates.canonicalize(Dates.CompoundPeriod(tend - tstart)))")
    return nothing
end
