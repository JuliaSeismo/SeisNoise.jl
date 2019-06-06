import Base:in, +, -, *, ==, convert, isempty, isequal, length, push!, sizeof, append!
import SeisIO: GeoLoc, PZResp
export CorrData

const corrfields = [:name, :loc, :comp, :rotated, :corr_type, :fs, :gain,
                    :freqmin, :freqmax, :cc_len, :cc_step, :whitened, :time_norm, :resp,
                    :maxlag]

"""
   CorrData

A structure for cross-correlations of ambient noise data.

## Fields: CorrData
| **Field** | **Description** |
|:------------|:-------------|
| :name       | Freeform channel names |
| :id         | Channel ids. use NET.STA.LOC.CHAN format when possible. |
| :loc        | Location (position) object. |
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
| :maxlag     | Maximum lag time in seconds to keep from correlations. |
| :t          | Starttime of each correlation. |
| :corr       | Correlations stored in columns. |
"""
mutable struct CorrData
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
  time_norm::Union{Bool,String}               # time normaliation
  resp::PZResp                                # response poles/zeros
  misc::Dict{String,Any}                      # misc
  notes::Array{String,1}                      # notes
  maxlag::Float64                             # maximum lag time [s]
  t::Array{Float64,1}                         # time
  corr::Array{<:Union{Float32,Float64},2}     # corr data

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
      time_norm::Union{Bool,String},
      resp     ::PZResp,
      misc     ::Dict{String,Any},
      notes    ::Array{String,1},
      maxlag   ::Float64,
      t        ::Array{Float64,1},
      corr     ::Array{<:Union{Float32,Float64},2}
      )

      return new(name, id, loc, comp, rotated, corr_type, fs, gain, freqmin,
                 freqmax, cc_len, cc_step, whitened, time_norm, resp, misc,
                 notes, maxlag, t, corr)
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
          time_norm::Union{Bool,String}        = false,
          resp     ::PZResp                    = PZResp(),
          misc     ::Dict{String,Any}          = Dict{String,Any}(),
          notes    ::Array{String,1}           = Array{String,1}(undef, 0),
          maxlag   ::Float64                   = zero(Float64),
          t        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          corr     ::Array{<:Union{Float32,Float64},2} = Array{Float32,2}(undef, 0, 2)
          ) = CorrData(name, id, loc, comp, rotated, corr_type, fs, gain,
                      freqmin, freqmax, cc_len, cc_step, whitened, time_norm,
                      resp, misc, notes, maxlag, t, corr)

CorrData(F1::FFTData, F2::FFTData, comp::String, rotated::Bool, corr_type::String,
        maxlag::Float64, t::Array{Float64,1}, corr::Array{<:Union{Float32,Float64},2}
       ) = CorrData(F1.name * '.' * F2.name, F1.id, F1.loc, comp, rotated,
                  corr_type, F1.fs, F1.gain, F1.freqmin, F1.freqmax, F1.cc_len,
                  F1.cc_step, F1.whitened, F1.time_norm, F1.resp, F1.misc,
                  F1.notes, maxlag, t, corr)

in(s::String, C::CorrData) = C.id==s

isempty(C::CorrData) = min(C.name=="",C.id == "",isempty(C.loc), C.comp=="",
                          C.rotated==false,C.fs == zero(typeof(C.fs)),
                          C.corr_type=="",
                          C.gain == one(typeof(C.gain)),
                          C.freqmin == zero(typeof(C.freqmin)),
                          C.freqmax == zero(typeof(C.freqmax)),
                          C.cc_len == zero(typeof(C.cc_len)),
                          C.cc_step == zero(typeof(C.cc_step)),
                          C.whitened == false,
                          C.time_norm == false,
                          isempty(C.resp), isempty(C.misc),isempty(C.notes),
                          C.maxlag == zero(typeof(C.maxlag)), isempty(C.t),
                          isempty(C.corr))

isequal(C::CorrData, U::CorrData) = minimum([hash(getfield(C,i))==hash(getfield(U,i)) for i in corrfields]::Array{Bool,1})
==(C::CorrData, U::CorrData) = isequal(C,U)::Bool

append!(C::CorrData, U::CorrData)  = (if C == U;
  setfield!(C, :t, vcat(getfield(C,:t), getfield(U,:t)));
  setfield!(C, :corr, hcat(getfield(C,:corr), getfield(U,:corr)))
  elseif isempty(C)
  [setfield!(C, i, getfield(U,i)) for i in corrfields];
  [setfield!(C, i, getfield(U,i)) for i in [:t,:corr]];
  end;
  return C)
+(C::CorrData, U::CorrData) = (T = deepcopy(C); return append!(T, U))

function sizeof(C::CorrData)
s = sum([sizeof(getfield(C,f)) for f in corrfields])
if !isempty(C.notes)
  s += sum([sizeof(i) for i in C.notes])
end
if !isempty(C.misc)
  s += sum([sizeof(i) for i in values(C.misc)])
end
return s
end
