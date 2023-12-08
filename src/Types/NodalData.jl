export gpu,cpu, NodalCorrData, NodalFFTData, NodalProcessedData


# functions for adapting NodalData to GPU
cpu(m::NodalData) = adapt(Array,m)
gpu(x::NodalData) = use_cuda[] ? cu(x) : x

# function for adapting NodalCorr to GPU
cpu(m::NodalCorrData) = adapt(Array,m)
gpu(x::NodalCorrData) = use_cuda[] ? cu(x) : x

Adapt.adapt_structure(to, N::NodalData) = NodalData(
    N.n,
    N.ox,
    N.oy,
    N.oz,
    N.info,
    N.id,
    N.name,
    N.loc,
    N.fs,
    N.gain,
    N.resp,
    N.units,
    N.src,
    N.misc,
    N.notes,
    N.t,
    adapt(to,N.data),
)

Adapt.adapt_structure(to, NC::NodalCorrData) = NodalCorrData(
    NC.n,
    NC.ox,
    NC.oy,
    NC.oz,
    NC.info,
    NC.id,
    NC.name,
    NC.loc,
    NC.fs,
    NC.gain,
    NC.freqmin,
    NC.freqmax,
    NC.cc_len,
    NC.maxlag,
    NC.time_norm,
    NC.whitened,
    NC.resp,
    NC.units,
    NC.src,
    NC.misc,
    NC.notes,
    NC.t,
    adapt(to,NC.corr),
)


mutable struct NodalCorrData <: NoiseData
    n::Int64
    ox::Float64                         # origin x
    oy::Float64                         # origin y
    oz::Float64                         # origin z
    info::Dict{String,Any}              # info
    id::Array{String,1}                 # id
    name::Array{String,1}               # name
    loc::Array{InstrumentPosition,1}    # loc
    fs::Array{Float64,1}                # fs
    gain::Array{Float64,1}              # gain
    freqmin::Float64                    # minumum frequency [Hz]
    freqmax::Float64                    # maximum frequency [Hz]
    cc_len::Int64                       # window_length [S]
    maxlag::Int64
    time_norm::String                   # time normaliation
    whitened::Bool
    resp::Array{InstrumentResponse,1}   # resp
    units::Array{String,1}              # units
    src::Array{String,1}                # src
    misc::Array{Dict{String,Any},1}     # misc
    notes::Array{Array{String,1},1}     # notes
    t::Array{Array{Int64,2},1}          # time
    corr::AbstractArray{<:AbstractFloat,2}  # fft data

    function NodalCorrData(
        n::Int64,
        ox::Float64,                         # origin x
        oy::Float64,                         # origin y
        oz::Float64,                         # origin z
        info::Dict{String,Any},              # info
        id::Array{String,1},                 # id
        name::Array{String,1},               # name
        loc::Array{InstrumentPosition,1},    # loc
        fs::Array{Float64,1},                # fs
        gain::Array{Float64,1},              # gain
        freqmin::Float64,                    # minumum frequency [Hz]
        freqmax::Float64,                    # maximum frequency [Hz]
        cc_len::Int64,                       # window_length [S]
        maxlag::Int64,                       # window_length [S]
        time_norm::String,                   # time normaliation
        whitened::Bool,
        resp::Array{InstrumentResponse,1},   # resp
        units::Array{String,1},              # units
        src::Array{String,1},                # src
        misc::Array{Dict{String,Any},1},     # misc
        notes::Array{Array{String,1},1},     # notes
        t::Array{Array{Int64,2},1},          # time
        corr::AbstractArray{<:AbstractFloat,2}  # fft data
        )

        return new(n,ox,oy,oz,info,id,name,loc,fs,gain,freqmin,freqmax,cc_len,
                   maxlag,time_norm,whitened,resp,units,src,misc,notes,t,corr)
    end

end


mutable struct NodalFFTData <: NoiseData
    n::Int64
    ns::Int64
    ox::Float64                         # origin x
    oy::Float64                         # origin y
    oz::Float64                         # origin z
    info::Dict{String,Any}              # info
    id::Array{String,1}                 # id
    name::Array{String,1}               # name
    loc::Array{InstrumentPosition,1}    # loc
    fs::Array{Float64,1}                # fs
    gain::Array{Float64,1}              # gain
    freqmin::Float64                    # minumum frequency [Hz]
    freqmax::Float64                    # maximum frequency [Hz]
    cc_len::Int64                       # window_length [S]
    time_norm::String                   # time normaliation
    resp::Array{InstrumentResponse,1}   # resp
    units::Array{String,1}              # units
    src::Array{String,1}                # src
    misc::Array{Dict{String,Any},1}     # misc
    notes::Array{Array{String,1},1}     # notes
    preprocessed::Bool
    dims::Array{Int64,1}                # FFT dimensions
    t::Array{Array{Int64,2},1}          # time
    fft::AbstractArray{<:Complex{<:AbstractFloat}}  # fft data

    function NodalFFTData(
        n::Int64,
        ns::Int64,
        ox::Float64,                         # origin x
        oy::Float64,                         # origin y
        oz::Float64,                         # origin z
        info::Dict{String,Any},              # info
        id::Array{String,1},                 # id
        name::Array{String,1},               # name
        loc::Array{InstrumentPosition,1},    # loc
        fs::Array{Float64,1},                # fs
        gain::Array{Float64,1},              # gain
        freqmin::Float64,                    # minumum frequency [Hz]
        freqmax::Float64,                    # maximum frequency [Hz]
        cc_len::Int64,                       # window_length [S]
        time_norm::String,                   # time normaliation
        resp::Array{InstrumentResponse,1},   # resp
        units::Array{String,1},              # units
        src::Array{String,1},                # src
        misc::Array{Dict{String,Any},1},     # misc
        notes::Array{Array{String,1},1},     # notes
        preprocessed::Bool,
        dims::Array{Int64,1},                # FFT dimensions
        t::Array{Array{Int64,2},1},          # time
        fft::AbstractArray{<:Complex{<:AbstractFloat}}  # fft data
        )

        return new(n,ns,ox,oy,oz,info,id,name,loc,fs,gain,freqmin,freqmax,cc_len,
                   time_norm,resp,units,src,misc,notes,preprocessed,dims,t,fft)
    end

end


mutable struct NodalProcessedData <: NoiseData
    n::Int64
    ns::Int64
    ox::Float64                         # origin x
    oy::Float64                         # origin y
    oz::Float64                         # origin z
    info::Dict{String,Any}              # info
    id::Array{String,1}                 # id
    name::Array{String,1}               # name
    loc::Array{InstrumentPosition,1}    # loc
    fs::Array{Float64,1}                # fs
    gain::Array{Float64,1}              # gain
    freqmin::Float64                    # minumum frequency [Hz]
    freqmax::Float64                    # maximum frequency [Hz]
    cc_len::Int64                       # window_length [S]
    time_norm::String                   # time normaliation
    resp::Array{InstrumentResponse,1}   # resp
    units::Array{String,1}              # units
    src::Array{String,1}                # src
    misc::Array{Dict{String,Any},1}     # misc
    notes::Array{Array{String,1},1}     # notes
    t::Array{Array{Int64,2},1}          # time
    data::AbstractArray{Float32, 3}     # fft data

    function NodalProcessedData(
        n::Int64,
        ns::Int64,
        ox::Float64,                         # origin x
        oy::Float64,                         # origin y
        oz::Float64,                         # origin z
        info::Dict{String,Any},              # info
        id::Array{String,1},                 # id
        name::Array{String,1},               # name
        loc::Array{InstrumentPosition,1},    # loc
        fs::Array{Float64,1},                # fs
        gain::Array{Float64,1},              # gain
        freqmin::Float64,                    # minumum frequency [Hz]
        freqmax::Float64,                    # maximum frequency [Hz]
        cc_len::Int64,                      # window_length [S]
        time_norm::String,                   # time normaliation
        resp::Array{InstrumentResponse,1},   # resp
        units::Array{String,1},              # units
        src::Array{String,1},                # src
        misc::Array{Dict{String,Any},1},     # misc
        notes::Array{Array{String,1},1},     # notes
        t::Array{Array{Int64,2},1},          # time
        data::AbstractArray{Float32, 3}      # fft data
        )

        return new(n,ns,ox,oy,oz,info,id,name,loc,fs,gain,
                   freqmin,freqmax,cc_len,time_norm,resp,
                   units,src,misc,notes,t,data)
    end

end
