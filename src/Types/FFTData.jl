import Base:in, +, -, *, ==, convert, isempty, isequal, length, push!, sizeof, append!
import SeisIO: GeoLoc, PZResp
export FFTData

const fftfields = [:name, :loc, :fs, :gain, :freqmin, :freqmax, :cc_step,
                   :cc_len, :whitened, :time_norm, :resp]

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
mutable struct FFTData
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
  time_norm::Union{Bool,String}               # time normaliation
  resp::PZResp                         # response poles/zeros
  misc::Dict{String,Any}                      # misc
  notes::Array{String,1}                      # notes
  t::Array{Float64,1}                         # time
  fft::Array{<:Union{Complex{Float32},Complex{Float64}},2}  # fft data

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
      time_norm::Union{Bool,String},
      resp     ::PZResp,
      misc     ::Dict{String,Any},
      notes    ::Array{String,1},
      t        ::Array{Float64,1},
      fft      ::Array{<:Union{Complex{Float32},Complex{Float64}},2}
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
          time_norm::Union{Bool,String}        = false,
          resp     ::PZResp                    = PZResp(),
          misc     ::Dict{String,Any}          = Dict{String,Any}(),
          notes    ::Array{String,1}           = Array{String,1}(undef, 0),
          t        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          fft      ::Array{<:Union{Complex{Float32},Complex{Float64}},2} =
                     Array{Complex{Float32},2}(undef, 0, 2)
          ) = FFTData(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
                     whitened, time_norm, resp, misc, notes, t, fft)

FFTData(C::SeisChannel,freqmin::Float64, freqmax::Float64,cc_len::Int,
        cc_step::Int, whitened::Bool,time_norm::Union{Bool,String},
        t::Array{Float64,1},fft::Array{<:Union{Complex{Float32},Complex{Float64}},2}
       ) = FFTData(C.id, Dates.format(u2d(C.t[1,2]*1e-6),"Y-mm-dd"), C.loc, C.fs, C.gain, freqmin, freqmax, cc_len,
                   cc_step, whitened, time_norm, C.resp, C.misc, C.notes, t, fft)

in(s::String, F::FFTData) = F.id==s

isempty(F::FFTData) = min(F.name=="",F.id == "",isempty(F.loc),
                          F.fs == zero(typeof(F.fs)),
                          F.gain == one(typeof(F.gain)),
                          F.freqmin == zero(typeof(F.freqmin)),
                          F.freqmax == zero(typeof(F.freqmax)),
                          F.cc_len == zero(typeof(F.cc_len)),
                          F.cc_step == zero(typeof(F.cc_step)),
                          F.whitened == false,
                          F.time_norm == false,
                          isempty(F.resp), isempty(F.misc),isempty(F.notes),
                          isempty(F.t),isempty(F.fft))

isequal(F::FFTData, U::FFTData) = minimum([hash(getfield(F,i))==hash(getfield(U,i)) for i in fftfields]::Array{Bool,1})
==(F::FFTData, U::FFTData) = isequal(F,U)::Bool

append!(F::FFTData, U::FFTData)  = (if F == U;
  setfield!(F, :t, vcat(getfield(F,:t), getfield(U,:t)));
  setfield!(F, :fft, hcat(getfield(F,:fft), getfield(U,:fft)))
  elseif isempty(F);
  [setfield!(F, i, getfield(U,i)) for i in fftfields];
  [setfield!(F, i, getfield(U,i)) for i in [:t,:fft]];
  end;
  return F)
+(F::FFTData, U::FFTData) = (T = deepcopy(F); return append!(T, U))

function sizeof(F::FFTData)
s = sum([sizeof(getfield(F,f)) for f in fftfields])
if !isempty(F.notes)
  s += sum([sizeof(i) for i in F.notes])
end
if !isempty(F.misc)
  s += sum([sizeof(i) for i in values(F.misc)])
end
return s
end
