import Base: in, +, -, *, ==, convert, isempty, isequal, length, push!, sizeof, append!
import Base: lastindex, firstindex, getindex, setindex!, sort, sort!
import SeisBase: GeoLoc, PZResp, InstrumentResponse
export RawData, FFTData, CorrData, gpu, cpu, sort, sort!, append!
export isequal, isempty, length, +, ==

# Many thanks to SeisBase for inspiration for NoiseData objects

const datafields = (:x,:fft,:corr)

abstract type NoiseData end

"""
  RawData

A structure for raw ambient noise data.

## Fields: RawData
| **Field** | **Description** |
|:--------------------|:--------------------|
| :name       | Freeform channel names |
| :id         | Channel ids. use NET.STA.LOC.CHAN format when possible. |
| :loc        | Location (position) object |
| :fs         | Sampling frequency in Hz. |
| :gain       | Scalar gain; divide data by the gain to convert to units  |
| :freqmin    | Minimum frequency from instrument response.  |
| :freqmax    | Maximum frequency from instrument response. |
| :cc_len     | Raw data window length in seconds. |
| :cc_step    | Spacing between windows in seconds. |
| :time_norm  | Apply one-bit whitening with "one_bit". |
| :resp       | SeisBase InstrumentResponse object |
| :misc       | Dictionary for non-critical information. |
| :notes      | Timestamped notes; includes automatically-logged acquisition and |
|             | processing information. |
| :t          | Starttime of each window. |
| :x          | Raw data stored in columns. |
"""
mutable struct RawData <: NoiseData
   name::String                                # name [Net.Sta.Loc.Chan]
   id::String                                  # id [Y-mm-dd] this is date of raw
   loc::GeoLoc                                 # loc
   fs::Float64                                 # sampling rate [Hz]
   gain::Float64                               # gain
   freqmin::Float64                            # minumum frequency [Hz]
   freqmax::Float64                            # maximum frequency [Hz]
   cc_len::Float64                             # window_length [S]
   cc_step::Float64                            # step between windows [S]
   whitened::Bool                              # whitening applied
   time_norm::String                           # time normaliation
   resp::InstrumentResponse                    # response poles/zeros
   misc::Dict{String,Any}                      # misc
   notes::Array{String,1}                      # notes
   t::Array{Float64,1}                         # time
   x::AbstractArray{<:AbstractFloat,2}           # raw data

 function RawData(
     name     ::String,
     id       ::String,
     loc      ::GeoLoc,
     fs       ::Float64,
     gain     ::Float64,
     freqmin  ::Float64,
     freqmax  ::Float64,
     cc_len   ::Float64,
     cc_step  ::Float64,
     whitened ::Bool,
     time_norm::String,
     resp     ::InstrumentResponse,
     misc     ::Dict{String,Any},
     notes    ::Array{String,1},
     t        ::Array{Float64,1},
     x        ::AbstractArray{<:AbstractFloat,2}
     )

     return new(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
             whitened, time_norm, resp, misc, notes, t, x)
   end

   function RawData(S::SeisData,cc_len::Real,cc_step::Real)
     merge!(S)
     # check cc_len & fs compatibility
     if mod(cc_len * S.fs[1], 1.) !== 0.
          throw(DomainError(cc_len, "cc_len * fs must be an integer for slicing."))
     end

     # check cc_step and sampling rate
     if mod(cc_step * S.fs[1], 1.) !== 0.
          throw(DomainError(cc_step, "cc_step * fs must be an integer for slicing."))
     end
     ungap!(S)

     # convert input to float 64
     cc_len = Float64(cc_len)
     cc_step = Float64(cc_step)

     # check if waveform length is < cc_len
     if length(S.x[1]) < cc_len * S.fs[1]
         throw(DomainError(cc_len, "cc_len must be less than length of data."))
     end

     # phase shift data onto exaclty the sampling rate
     phase_shift!(S)
     # subset by time
     startslice, endslice = nearest_start_end(S[1],cc_len, cc_step)
     # get starting and ending indices
     startind, endind = slide_ind(startslice, endslice, S.fs[1], S.t[1])
     x, starts = slide(@view(S.x[1][startind:endind]), cc_len, cc_step, S.fs[1], startslice,endslice)
     return new(S[1].id,Dates.format(u2d(starts[1]),"Y-mm-dd"),S.loc[1],
                S.fs[1],S.gain[1],1. / cc_len,S[1].fs/2,cc_len,cc_step,false,
                "",S[1].resp,S.misc[1],S.notes[1],starts,x)
   end

   function RawData(C::SeisChannel,cc_len::Real,cc_step::Real)
    # check cc_len and sampling rate
    if mod(cc_len * C.fs, 1.) !== 0.
         throw(DomainError(cc_len, "cc_len * fs must be an integer for slicing."))
    end

    # check cc_step and sampling rate
    if mod(cc_step * C.fs, 1.) !== 0.
         throw(DomainError(cc_step, "cc_step * fs must be an integer for slicing."))
    end

    ungap!(C)
    # check if waveform length is < cc_len
    if length(C.x) < cc_len * C.fs
        throw(DomainError(cc_len, "cc_len must be less than length of data."))
    end

    # convert input to float 64
    cc_len = Float64(cc_len)
    cc_step = Float64(cc_step)

    # phase shift data onto exaclty the sampling rate
    phase_shift!(C)

    # subset by time
    startslice, endslice = nearest_start_end(C,cc_len, cc_step)

    # get starting and ending indices
    startind, endind = slide_ind(startslice, endslice, C.fs,C.t)
    x, starts = slide(@view(C.x[startind:endind]), cc_len, cc_step, C.fs, startslice,endslice)
    return new(C.id,Dates.format(u2d(starts[1]),"Y-mm-dd"),C.loc,C.fs,
               C.gain,1. / cc_len,C.fs/2,cc_len,cc_step,false,"", C.resp,
               C.misc,C.notes,starts,x)
   end
end

RawData(;
         name     ::String                    = "",
         id       ::String                    = "",
         loc      ::GeoLoc                    = GeoLoc(),
         fs       ::Float64                   = zero(Float64),
         gain     ::Float64                   = one(Float64),
         freqmin  ::Float64                   = zero(Float64),
         freqmax  ::Float64                   = zero(Float64),
         cc_len   ::Float64                   = zero(Float64),
         cc_step  ::Float64                   = zero(Float64),
         whitened ::Bool                      = false,
         time_norm::String                    = "",
         resp     ::InstrumentResponse        = PZResp(),
         misc     ::Dict{String,Any}          = Dict{String,Any}(),
         notes    ::Array{String,1}           = Array{String,1}(undef, 0),
         t        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
         x        ::AbstractArray{<:AbstractFloat,2} = Array{Float32,2}(undef, 0, 2)
         ) = RawData(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
                whitened, time_norm, resp, misc, notes, t, x)

"""
FFTData

A structure for fourier transforms (FFT) of ambient noise data.

## Fields: FFTData
| **Field** | **Description** |
|:--------------------|:--------------------|
| :name       | Freeform channel names |
| :id         | Channel ids. use NET.STA.LOC.CHAN format when possible. |
| :loc        | Location (position) object |
| :fs         | Sampling frequency in Hz. |
| :gain       | Scalar gain; divide data by the gain to convert to units  |
| :freqmin    | Minimum frequency for whitening.  |
| :freqmax    | Maximum frequency for whitening. |
| :cc_len     | Raw data window length in number of points. |
| :cc_step    | Spacing between windows in number of points. |
| :whitened   | Whitening applied.
| :time_norm  | Apply one-bit whitening with "one_bit". |
| :resp       | SeisBase InstrumentResponse object |
| :misc       | Dictionary for non-critical information. |
| :notes      | Timestamped notes; includes automatically-logged acquisition and |
|             | processing information. |
| :t          | Starttime of each FFT. |
| :fft        | FFTs stored in columns. |
"""
mutable struct FFTData <: NoiseData
   name::String                                # name [Net.Sta.Loc.Chan]
   id::String                                  # id [Y-mm-dd] this is date of fft
   loc::GeoLoc                                 # loc
   fs::Float64                                 # sampling rate [Hz]
   gain::Float64                               # gain
   freqmin::Float64                            # minumum frequency [Hz]
   freqmax::Float64                            # maximum frequency [Hz]
   cc_len::Float64                             # window_length [S]
   cc_step::Float64                            # step between windows [S]
   whitened::Bool                              # whitening applied
   time_norm::String                           # time normaliation
   resp::InstrumentResponse                    # response poles/zeros
   misc::Dict{String,Any}                      # misc
   notes::Array{String,1}                      # notes
   t::Array{Float64,1}                         # time
   fft::AbstractArray{<:Complex{<:AbstractFloat},2}  # fft data

function FFTData(
   name     ::String,
   id       ::String,
   loc      ::GeoLoc,
   fs       ::Float64,
   gain     ::Float64,
   freqmin  ::Float64,
   freqmax  ::Float64,
   cc_len   ::Float64,
   cc_step  ::Float64,
   whitened ::Bool,
   time_norm::String,
   resp     ::InstrumentResponse,
   misc     ::Dict{String,Any},
   notes    ::Array{String,1},
   t        ::Array{Float64,1},
   fft      ::AbstractArray{<:Complex{<:AbstractFloat},2}
   )

   return new(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
           whitened, time_norm, resp, misc, notes, t, fft)
 end
end

FFTData(;
       name     ::String                    = "",
       id       ::String                    = "",
       loc      ::GeoLoc                    = GeoLoc(),
       fs       ::Float64                   = zero(Float64),
       gain     ::Float64                   = one(Float64),
       freqmin  ::Float64                   = zero(Float64),
       freqmax  ::Float64                   = zero(Float64),
       cc_len   ::Float64                   = zero(Float64),
       cc_step  ::Float64                   = zero(Float64),
       whitened ::Bool                      = false,
       time_norm::String                    = "",
       resp     ::InstrumentResponse        = PZResp(),
       misc     ::Dict{String,Any}          = Dict{String,Any}(),
       notes    ::Array{String,1}           = Array{String,1}(undef, 0),
       t        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
       fft      ::AbstractArray{<:Complex{<:AbstractFloat},2} =
                  Array{Complex{Float32},2}(undef, 0, 2)
       ) = FFTData(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
             whitened, time_norm, resp, misc, notes, t, fft)

# FFTData(C::SeisChannel,freqmin::Float64, freqmax::Float64,cc_len::Int,
#      cc_step::Int, whitened::Bool,time_norm::Union{Bool,String},
#      t::Array{Float64,1},fft::Array{<:Union{Complex{Float32},Complex{Float64}},2}
#     ) = FFTData(C.id, Dates.format(u2d(C.t[1,2]*1e-6),"Y-mm-dd"), C.loc, C.fs, C.gain, freqmin, freqmax, cc_len,
#                 cc_step, whitened, time_norm, C.resp, C.misc, C.notes, t, fft)

"""
   CorrData

A structure for cross-correlations of ambient noise data.

## Fields: CorrData
| **Field** | **Description** |
|:------------|:-------------|
| :name       | Freeform channel name in NET1.STA1.LOC1.CHAN1.NET2.STA2.LOC2.CHAN2 form |
| :id         | Channel ids. use NET.STA.LOC.CHAN format when possible. |
| :loc        | Location (position) object for station 1. |
| :comp       | Component of cross-correlation (e.g ZZ, RT, etc..). |
| :rotated    | Boolean for rotation. |
| :fs         | Sampling frequency in Hz. |
| :gain       | Scalar gain; divide data by the gain to convert to units.  |
| :freqmin    | Minimum frequency for whitening.  |
| :freqmax    | Maximum frequency for whitening. |
| :cc_len     | Raw data window length in number of points. |
| :cc_step    | Spacing between windows in number of points. |
| :whitened   | Whitening applied.
| :time_norm  | Apply one-bit whitening with "one_bit". |
| :resp       | SeisBase InstrumentResponse object |
| :misc       | Dictionary for non-critical information. |
| :notes      | Timestamped notes; includes automatically-logged acquisition and |
|             | processing information. |
| :dist       | Distance between station 1 and station 2 in Km. |
| :azi        | Azimuth from station 1 and station 2 in degrees. |
| :baz        | Back azimuth between station 1 and station 2 in degrees. |
| :maxlag     | Maximum lag time in seconds to keep from correlations. |
| :t          | Starttime of each correlation. |
| :corr       | Correlations stored in columns. |
"""
mutable struct CorrData <: NoiseData
  name::String                                # name [Net1.Sta1.Loc1.Chan1.Net2.Sta2.Loc2.Chan2]
  id::String                                  # id [Y-mm-dd] this is date of corr
  loc::GeoLoc                                 # loc
  comp::String                                # 1st channel, 2nd channel [ZZ,RT,..]
  rotated::Bool                               # rotation applied?
  corr_type::String                           # "corr", "deconv", "coh"
  fs::Float64                                 # sampling rate [Hz]
  gain::Float64                               # gain
  freqmin::Float64                            # minumum frequency [Hz]
  freqmax::Float64                            # maximum frequency [Hz]
  cc_len::Float64                             # window_length [S]
  cc_step::Float64                            # step between windows [S]
  whitened::Bool                              # whitening applied
  time_norm::String                           # time normaliation
  resp::InstrumentResponse                    # response poles/zeros
  misc::Dict{String,Any}                      # misc
  notes::Array{String,1}                      # notes
  dist::Float64                               # distance between stations [km]
  azi::Float64                                # azimuth between stations [degrees]
  baz::Float64                                # back azimuth between stations [degrees]
  maxlag::Float64                             # maximum lag time [s]
  t::Array{Float64,1}                         # time
  corr::AbstractArray{<:AbstractFloat,2}     # corr data

  function CorrData(
      name     ::String,
      id       ::String,
      loc      ::GeoLoc,
      comp     ::String,
      rotated  ::Bool,
      corr_type::String,
      fs       ::Float64,
      gain     ::Float64,
      freqmin  ::Float64,
      freqmax  ::Float64,
      cc_len   ::Float64,
      cc_step  ::Float64,
      whitened ::Bool,
      time_norm::String,
      resp     ::InstrumentResponse,
      misc     ::Dict{String,Any},
      notes    ::Array{String,1},
      dist     ::Float64,
      azi      ::Float64,
      baz      ::Float64,
      maxlag   ::Float64,
      t        ::Array{Float64,1},
      corr     ::AbstractArray{<:AbstractFloat,2}
      )

      return new(name, id, loc, comp, rotated, corr_type, fs, gain, freqmin,
              freqmax, cc_len, cc_step, whitened, time_norm, resp, misc,
              notes, dist, azi, baz, maxlag, t, corr)
    end
end

CorrData(;
          name     ::String                    = "",
          id       ::String                    = "",
          loc      ::GeoLoc                    = GeoLoc(),
          comp     ::String                    = "",
          rotated  ::Bool                      = false,
          corr_type::String                    = "",
          fs       ::Float64                   = zero(Float64),
          gain     ::Float64                   = one(Float64),
          freqmin  ::Float64                   = zero(Float64),
          freqmax  ::Float64                   = zero(Float64),
          cc_len   ::Float64                   = zero(Float64),
          cc_step  ::Float64                   = zero(Float64),
          whitened ::Bool                      = false,
          time_norm::String                    = "",
          resp     ::InstrumentResponse        = PZResp(),
          misc     ::Dict{String,Any}          = Dict{String,Any}(),
          notes    ::Array{String,1}           = Array{String,1}(undef, 0),
          dist     ::Float64                   = zero(Float64),
          azi      ::Float64                   = zero(Float64),
          baz      ::Float64                   = zero(Float64),
          maxlag   ::Float64                   = zero(Float64),
          t        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          corr     ::AbstractArray{<:AbstractFloat,2} = Array{Float32,2}(undef, 0, 2)
          ) = CorrData(name, id, loc, comp, rotated, corr_type, fs, gain,
                 freqmin, freqmax, cc_len, cc_step, whitened, time_norm,
                 resp, misc, notes, dist, azi, baz, maxlag, t, corr)

CorrData(F1::FFTData, F2::FFTData, comp::String, rotated::Bool, corr_type::String,
       maxlag::Float64, t::Array{Float64,1}, corr::AbstractArray{<:AbstractFloat,2}
      ) = CorrData(F1.name * '.' * F2.name, F1.id, F1.loc, comp, rotated,
                 corr_type, F1.fs, F1.gain, F1.freqmin, F1.freqmax, F1.cc_len,
                 F1.cc_step, F1.whitened, F1.time_norm, F1.resp, F1.misc,
                 F1.notes, get_dist(F1,F2), get_azi(F1,F2), get_baz(F1,F2),
                 maxlag, t, corr)


in(s::String, A::NoiseData) = A.id==s
length(N::NoiseData) = length(getfield(N,:t))
firstindex(N::NoiseData) = 1
lastindex(N::NoiseData) = length(N)

function getindex(N::NoiseData, J::Array{Int,1})
  n = length(N)
  T = typeof(N)
  U = T()
  F = fieldnames(T)
  for f in F
     if f in datafields
        setfield!(U, f, getindex(getfield(N,f),:,J))
     elseif f == :t
        setfield!(U, f, getindex(getfield(N,f),J))
     else
        setfield!(U, f, getfield(N,f))
     end
  end
  return U
end

function getindex(N::NoiseData, J::Array{DateTime,1})
   n = length(N)
   T = typeof(N)
   U = T()
   F = fieldnames(T)
   t = d2u.(J)
   ind = findall(in(t),N.t)
   for f in F
      if f in datafields
         setfield!(U, f, getindex(getfield(N,f),:,ind))
      elseif f == :t
         setfield!(U, f, getindex(getfield(N,f),ind))
      else
         setfield!(U, f, getfield(N,f))
      end
   end
   return U
 end
getindex(N::NoiseData, J::AbstractRange) = getindex(N, collect(J))

function setindex!(N::T, U::T, J::Array{Int,1}) where {T<:NoiseData}
  typeof(N) == typeof(U) || throw(MethodError)
  length(J) == length(U) || throw(BoundsError)
  length(N) >= maximum(J) || throw(BoundsError)
  F = fieldnames(T)
  for f in F
     if f in datafields
        setindex!(getfield(N, f), getfield(U, f), :,J)
     elseif f == :t
        setindex!(getfield(N,f),getfield(U,f),J)
     end
  end
  # ([(getfield(S, f))[J] = getfield(U, f) for f in datafields])
  return nothing
end
setindex!(N::NoiseData, U::NoiseData, J::AbstractRange) = setindex!(N, U, collect(J))

isempty(R::RawData) = min(R.name=="",R.id == "",isempty(R.loc),
                          R.fs == zero(typeof(R.fs)),
                          R.gain == one(typeof(R.gain)),
                          R.freqmin == zero(typeof(R.freqmin)),
                          R.freqmax == zero(typeof(R.freqmax)),
                          R.cc_len == zero(typeof(R.cc_len)),
                          R.cc_step == zero(typeof(R.cc_step)),
                          R.time_norm == "",
                          isempty(R.resp), isempty(R.misc),isempty(R.notes),
                          isempty(R.t),isempty(R.x))

isempty(F::FFTData) = min(F.name=="",F.id == "",isempty(F.loc),
                         F.fs == zero(typeof(F.fs)),
                         F.gain == one(typeof(F.gain)),
                         F.freqmin == zero(typeof(F.freqmin)),
                         F.freqmax == zero(typeof(F.freqmax)),
                         F.cc_len == zero(typeof(F.cc_len)),
                         F.cc_step == zero(typeof(F.cc_step)),
                         F.whitened == false,
                         F.time_norm == "",
                         isempty(F.resp), isempty(F.misc),isempty(F.notes),
                         isempty(F.t),isempty(F.fft))

isempty(C::CorrData) = min(C.name=="",C.id == "",isempty(C.loc), C.comp=="",
                        C.rotated==false,C.fs == zero(typeof(C.fs)),
                        C.corr_type=="",
                        C.gain == one(typeof(C.gain)),
                        C.freqmin == zero(typeof(C.freqmin)),
                        C.freqmax == zero(typeof(C.freqmax)),
                        C.cc_len == zero(typeof(C.cc_len)),
                        C.cc_step == zero(typeof(C.cc_step)),
                        C.whitened == false,
                        C.time_norm == "",
                        isempty(C.resp), isempty(C.misc),isempty(C.notes),
                        C.maxlag == zero(typeof(C.maxlag)), isempty(C.t),
                        isempty(C.corr))

function isequal(N::T,U::T) where {T<:NoiseData}
   return minimum([hash(getfield(N,i))==hash(getfield(U,i)) for i in fieldnames(T)]::Array{Bool,1})
end
==(N::T,U::T) where {T<:NoiseData} = isequal(N,U)

function append!(N::T, U::T) where {T<:NoiseData}
   F = fieldnames(T)
   if N.name == U.name
      d = intersect(datafields,F)[1]
      setfield!(N, :t, vcat(getfield(N,:t), getfield(U,:t)))
      setfield!(N, d, hcat(getfield(N,d), getfield(U,d)))
   elseif isempty(N)
      [setfield!(N, f, getfield(U,f)) for f in F]
   end
   return N
end

function sort!(N::T; rev=false::Bool) where {T<:NoiseData}
  j = sortperm(getfield(N, :t), rev=rev)
  F = fieldnames(T)
  d = intersect(datafields,F)[1]
  setfield!(N, d, getfield(N,d)[:,j])
  setfield!(N, :t, getfield(N,:t)[j])
  return nothing
end

function sort(N::T; rev=false::Bool) where {T<:NoiseData}
  U = deepcopy(N)
  sort!(U, rev=rev)
  return U
end

+(N::T, U::T) where {T<:NoiseData} = (Z = deepcopy(N); return append!(Z, U))

# Allow for transfer from CPU to GPU

Adapt.adapt_structure(to, R::RawData) = RawData(R.name, R.id, R.loc, R.fs, R.gain,
                                                R.freqmin, R.freqmax, R.cc_len,
                                                R.cc_step, R.whitened, R.time_norm,
                                                R.resp, R.misc, R.notes, R.t,
                                                adapt(to, R.x))
Adapt.adapt_structure(to, F::FFTData) = FFTData(F.name, F.id, F.loc, F.fs, F.gain,
                                                F.freqmin, F.freqmax, F.cc_len,
                                                F.cc_step, F.whitened, F.time_norm,
                                                F.resp, F.misc, F.notes, F.t,
                                                adapt(to, F.fft))
Adapt.adapt_structure(to, C::CorrData) = CorrData(C.name, C.id, C.loc, C.comp,
                                                  C.rotated, C.corr_type, C.fs,
                                                  C.gain, C.freqmin, C.freqmax,
                                                  C.cc_len, C.cc_step, C.whitened,
                                                  C.time_norm, C.resp, C.misc,
                                                  C.notes, C.dist,C.azi,C.baz,
                                                  C.maxlag,C.t,
                                                  adapt(to, C.corr))

# CPU/GPU movement conveniences
cpu(m::NoiseData) = adapt(Array,m)
gpu(x::NoiseData) = use_cuda[] ? cu(x) : x
