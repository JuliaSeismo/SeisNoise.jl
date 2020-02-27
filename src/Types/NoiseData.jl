import Base:in, +, -, *, ==, convert, isempty, isequal, length, push!, sizeof, append!
import SeisIO: GeoLoc, PZResp
export RawData, FFTData, CorrData, gpu, cpu

const rawfields = [:name, :loc, :fs, :gain, :freqmin, :freqmax, :cc_step,
                   :cc_len, :whitened, :time_norm, :resp, :x]
const fftfields = [:name, :loc, :fs, :gain, :freqmin, :freqmax, :cc_step,
                   :cc_len, :whitened, :time_norm, :resp, :fft]

const corrfields = [:name, :loc, :comp, :rotated, :corr_type, :fs, :gain,
                   :freqmin, :freqmax, :cc_len, :cc_step, :whitened, :time_norm, :resp,
                   :maxlag, :dist, :azi, :baz, :corr]

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
| :resp       | Instrument response object, format [zeros poles] |
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
   cc_len::Int                                 # window_length [s]
   cc_step::Int                                # step between windows [s]
   whitened::Bool                              # whitening applied
   time_norm::String                           # time normaliation
   resp::PZResp                                # response poles/zeros
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
     cc_len   ::Int,
     cc_step  ::Int,
     whitened ::Bool,
     time_norm::String,
     resp     ::PZResp,
     misc     ::Dict{String,Any},
     notes    ::Array{String,1},
     t        ::Array{Float64,1},
     x        ::AbstractArray{<:AbstractFloat,2}
     )

     return new(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
             whitened, time_norm, resp, misc, notes, t, x)
   end

   function RawData(S::SeisData,cc_len::Int,cc_step::Int)
     merge!(S)
     # subset by time
     startslice, endslice = u2d.(nearest_start_end(S[1],cc_len, cc_step))
       # check if waveform length is < cc_len
     if Int(floor((endslice - startslice).value / 1000)) < cc_len
         return nothing
     end

     # phase shift data onto exaclty the sampling rate
     phase_shift!(S)
     # get starting and ending indices
     startind, endind = slide_ind(startslice, endslice, S[1].t)
     x, starts = slide(S[1].x[startind:endind] , cc_len, cc_step, S[1].fs, S[1].t)
     return new(S[1].id,Dates.format(u2d(S[1].t[1,2]*1e-6),"Y-mm-dd"),S[1].loc,
                S[1].fs,S[1].gain,1/cc_len,S[1].fs/2,cc_len,cc_step,false,
                "",S[1].resp,S[1].misc,S[1].notes,starts,x)
   end

   function RawData(C::SeisChannel,cc_len::Int,cc_step::Int)
     # subset by time
     startslice, endslice = u2d.(nearest_start_end(C,cc_len, cc_step))
       # check if waveform length is < cc_len
     if Int(floor((endslice - startslice).value / 1000)) < cc_len
         return nothing
     end

     # phase shift data onto exaclty the sampling rate
     phase_shift!(C)
     # get starting and ending indices
     startind, endind = slide_ind(startslice, endslice, C.t)
     x, starts = slide(C.x[startind:endind] , cc_len, cc_step, C.fs, C.t)
     return new(C.id,Dates.format(u2d(C.t[1,2]*1e-6),"Y-mm-dd"),C.loc,C.fs,
                C.gain,1/cc_len,C.fs/2,cc_len,cc_step,false,"", C.resp,
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
         cc_len   ::Int                       = zero(Int),
         cc_step  ::Int                       = zero(Int),
         whitened ::Bool                      = false,
         time_norm::String                    = "",
         resp     ::PZResp                    = PZResp(),
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
| :cc_len     | Length of each correlation in seconds. |
| :cc_step    | Spacing between correlation windows in seconds. |
| :whitened   | Whitening applied.
| :time_norm  | Apply one-bit whitening with "one_bit". |
| :resp       | Instrument response object, format [zeros poles] |
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
   cc_len::Int                                 # window_length [s]
   cc_step::Int                                # step between windows [s]
   whitened::Bool                              # whitening applied
   time_norm::String                           # time normaliation
   resp::PZResp                                # response poles/zeros
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
   cc_len   ::Int,
   cc_step  ::Int,
   whitened ::Bool,
   time_norm::String,
   resp     ::PZResp,
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
       cc_len   ::Int                       = zero(Int),
       cc_step  ::Int                       = zero(Int),
       whitened ::Bool                      = false,
       time_norm::String                    = "",
       resp     ::PZResp                    = PZResp(),
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
| :cc_len     | Length of each correlation in seconds. |
| :cc_step    | Spacing between correlation windows in seconds. |
| :whitened   | Whitening applied.
| :time_norm  | Apply one-bit whitening with "one_bit". |
| :resp       | Instrument response object, format [zeros poles] |
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
  cc_len::Int                                 # window_length [s]
  cc_step::Int                                # step between windows [s]
  whitened::Bool                              # whitening applied
  time_norm::String                           # time normaliation
  resp::PZResp                                # response poles/zeros
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
      cc_len   ::Int,
      cc_step  ::Int,
      whitened ::Bool,
      time_norm::String,
      resp     ::PZResp,
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
          cc_len   ::Int                       = zero(Int),
          cc_step  ::Int                       = zero(Int),
          whitened ::Bool                      = false,
          time_norm::String                    = "",
          resp     ::PZResp                    = PZResp(),
          misc     ::Dict{String,Any}          = Dict{String,Any}(),
          notes    ::Array{String,1}           = Array{String,1}(undef, 0),
          dist     ::Float64                   = zero(Float64),
          azi      ::Float64                   = zero(Float64),
          baz      ::Float64                   = zero(Float64),
          maxlag   ::Float64                   = zero(Float64),
          t        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          corr     ::AbstractArray{AbstractFloat,2} = Array{Float32,2}(undef, 0, 2)
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

isequal(R::RawData, U::RawData) = minimum([hash(getfield(R,i))==hash(getfield(U,i)) for i in rawfields]::Array{Bool,1})
==(R::RawData, U::RawData) = isequal(R,U)::Bool

isequal(F::FFTData, U::FFTData) = minimum([hash(getfield(F,i))==hash(getfield(U,i)) for i in fftfields]::Array{Bool,1})
==(F::FFTData, U::FFTData) = isequal(F,U)::Bool

isequal(C::CorrData, U::CorrData) = minimum([hash(getfield(C,i))==hash(getfield(U,i)) for i in corrfields]::Array{Bool,1})
==(C::CorrData, U::CorrData) = isequal(C,U)::Bool

append!(R::RawData, U::RawData)  = (if R.name == U.name;
   setfield!(R, :t, vcat(getfield(R,:t), getfield(U,:t)));
   setfield!(R, :x, hcat(getfield(R,:x), getfield(U,:x)))
   elseif isempty(R);
   [setfield!(R, i, getfield(U,i)) for i in rawfields];
   [setfield!(R, i, getfield(U,i)) for i in [:t,:x]];
   end;
   return R)

append!(F::FFTData, U::FFTData)  = (if F.name == U.name;
   setfield!(F, :t, vcat(getfield(F,:t), getfield(U,:t)));
   setfield!(F, :fft, hcat(getfield(F,:fft), getfield(U,:fft)))
   elseif isempty(F);
   [setfield!(F, i, getfield(U,i)) for i in fftfields];
   [setfield!(F, i, getfield(U,i)) for i in [:t,:fft]];
   end;
   return F)

 append!(C::CorrData, U::CorrData)  = (if C == U;
   setfield!(C, :t, vcat(getfield(C,:t), getfield(U,:t)));
   setfield!(C, :corr, hcat(getfield(C,:corr), getfield(U,:corr)))
   elseif isempty(C)
   [setfield!(C, i, getfield(U,i)) for i in corrfields];
   [setfield!(C, i, getfield(U,i)) for i in [:t,:corr]];
   end;
   return C)

+(R::RawData, U::RawData) = (T = deepcopy(R); return append!(T, U))
+(F::FFTData, U::FFTData) = (T = deepcopy(F); return append!(T, U))
+(C::CorrData, U::CorrData) = (T = deepcopy(C); return append!(T, U))

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
