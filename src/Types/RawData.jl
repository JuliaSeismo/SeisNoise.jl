import Base:in, +, -, *, ==, convert, isempty, isequal, length, push!, sizeof, append!
import SeisIO: GeoLoc, PZResp
export RawData

const rawfields = [:name, :loc, :fs, :gain, :freqmin, :freqmax, :cc_step,
                   :cc_len, :time_norm, :resp]

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
mutable struct RawData
  name::String                                # name [Net.Sta.Loc.Chan]
  id::String                                  # id [Y-mm-dd] this is date of raw
  loc::GeoLoc                                 # loc
  fs::Float64                                 # sampling rate [Hz]
  gain::Float64                               # gain
  freqmin::Float64                            # minumum frequency [Hz]
  freqmax::Float64                            # maximum frequency [Hz]
  cc_len::Int                                 # window_length [s]
  cc_step::Int                                # step between windows [s]
  time_norm::Union{Bool,String}               # time normaliation
  resp::PZResp                                # response poles/zeros
  misc::Dict{String,Any}                      # misc
  notes::Array{String,1}                      # notes
  t::Array{Float64,1}                         # time
  x::Array{<:Union{Float32,Float64},2}  # raw data

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
      time_norm::Union{Bool,String},
      resp     ::PZResp,
      misc     ::Dict{String,Any},
      notes    ::Array{String,1},
      t        ::Array{Float64,1},
      x        ::Array{<:Union{Float32,Float64},2}
      )

      return new(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
                time_norm, resp, misc, notes, t, x)
    end

    function RawData(S::SeisData,cc_len::Int,cc_step::Int)
      merge!(S)
      # subset by time
      starttime, endtime = u2d.(nearest_start_end(S[1],cc_len, cc_step))
        # check if waveform length is < cc_len
      if Int(floor((endtime - starttime).value / 1000)) < cc_len
          return nothing
      end
      sync!(S,s=starttime,t=endtime) # sync start and end times
      x, starts = slide(S[1], cc_len, cc_step)
      return new(S[1].name,S[1].id,S[1].loc,S[1].fs,S[1].gain,1/cc_len,S[1].fs/2,
                cc_len,cc_step,false,S[1].resp,S[1].misc,S[1].notes,starts,x)
    end

    function RawData(C::SeisChannel,cc_len::Int,cc_step::Int)
      # subset by time
      starttime, endtime = u2d.(nearest_start_end(C,cc_len, cc_step))
        # check if waveform length is < cc_len
      if Int(floor((endtime - starttime).value / 1000)) < cc_len
          return nothing
      end
      sync!(C,s=starttime,t=endtime) # sync start and end times
      x, starts = slide(C, cc_len, cc_step)
      return new(C.name,C.id,C.loc,C.fs,C.gain,1/cc_len,C.fs/2,
                cc_len,cc_step,false,C.resp,C.misc,C.notes,starts,x)
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
          time_norm::Union{Bool,String}        = false,
          resp     ::PZResp                    = PZResp(),
          misc     ::Dict{String,Any}          = Dict{String,Any}(),
          notes    ::Array{String,1}           = Array{String,1}(undef, 0),
          t        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          x      ::Array{<:Union{Float32,Float64},2} =
                     Array{Float32,2}(undef, 0, 2)
          ) = RawData(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
                     time_norm, resp, misc, notes, t, x)

in(s::String, R::RawData) = R.id==s

isempty(R::RawData) = min(R.name=="",R.id == "",isempty(R.loc),
                          R.fs == zero(typeof(F.fs)),
                          R.gain == one(typeof(F.gain)),
                          R.freqmin == zero(typeof(F.freqmin)),
                          R.freqmax == zero(typeof(F.freqmax)),
                          R.cc_len == zero(typeof(F.cc_len)),
                          R.cc_step == zero(typeof(F.cc_step)),
                          F.time_norm == false,
                          isempty(R.resp), isempty(R.misc),isempty(R.notes),
                          isempty(R.t),isempty(R.raw))

isequal(R::RawData, U::RawData) = minimum([hash(getfield(R,i))==hash(getfield(U,i)) for i in rawfields]::Array{Bool,1})
==(F::RawData, U::RawData) = isequal(R,U)::Bool

append!(R::RawData, U::RawData)  = (if R == U;
  setfield!(R, :t, vcat(getfield(R,:t), getfield(U,:t)));
  setfield!(R, :raw, hcat(getfield(R,:raw), getfield(U,:raw)))
  elseif isempty(R);
  [setfield!(R, i, getfield(U,i)) for i in rawfields];
  [setfield!(R, i, getfield(U,i)) for i in [:t,:raw]];
  end;
  return R)
+(R::RawData, U::RawData) = (T = deepcopy(R); return append!(T, U))

function sizeof(R::RawData)
s = sum([sizeof(getfield(R,f)) for f in rawfields])
if !isempty(R.notes)
  s += sum([sizeof(i) for i in R.notes])
end
if !isempty(R.misc)
  s += sum([sizeof(i) for i in values(R.misc)])
end
return s
end
