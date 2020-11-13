export gpu,cpu, NodalCorr

mutable struct NodalCorr <: NoiseData
    n::Int64
    id::Array{String,1}                 # id
    name::Array{String,1}               # name
    loc::Array{<:InstrumentPosition,1}  # loc
    fs::Float64                         # fs
    src::String                         # src
    t::Array{Array{Int64,2},1}          # time
    corr::AbstractArray{Float32, 2}     # actual data

    function NodalCorr(
        n        ::Int64,
        id       ::Array{String,1},
        name     ::Array{String,1},
        loc      ::Array{<:InstrumentPosition,1},
        fs       ::Float64,
        src      ::String,
        t        ::Array{Array{Int64,2},1} ,
        corr     ::AbstractArray{<:AbstractFloat,2},
        )
  
        return new(n, id, name, loc, fs, src, t, corr)
    end
end

NodalCorr(;
          n        ::Int64                       = 0,
          id       ::Array{String,1}             = Array{String}(undef,0),
          name     ::Array{String,1}             = Array{String}(undef,0),
          loc      ::Array{InstrumentPosition,1} = Array{InstrumentPosition}(undef,0), 
          fs       ::Float64                     = zero(Float64),
          src      ::String                      = "",
          t        ::Array{Array{Int64,2},1}     = Array{Array{Int64,2},1}(undef,0),
          corr     ::AbstractArray{<:AbstractFloat,2} = Array{Float32,2}(undef, 0, 2),
          ) = CorrData(n, id, name, loc, fs, src, t, corr)


# functions for adapting NodalData to GPU 
cpu(m::NodalData) = adapt(Array,m)
gpu(x::NodalData) = use_cuda[] ? cu(x) : x

# function for adapting NodalCorr to GPU 
cpu(m::NodalCorr) = adapt(Array,m)
gpu(x::NodalCorr) = use_cuda[] ? cu(x) : x

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

Adapt.adapt_structure(to, N::NodalCorr) = NodalCorr(
    N.n,
    N.id,
    N.name,
    N.loc,
    N.fs,
    N.gain,
    N.resp,
    N.units,
    N.src,
    N.t,
    adapt(to,N.corr),
)
