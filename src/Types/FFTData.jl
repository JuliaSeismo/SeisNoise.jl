import Base:in, +, -, *, ==, convert, isempty, isequal, length, push!, sizeof
export FFTData

const fftfields = [:id, :name, :loc, :fs, :gain, :freqmin, :freqmax, :cc_step,
                   :whitened, :time_norm, :resp,:notes, :misc, :t, :fft]

# This is type-stable for F = FFTData() but not for keyword args
mutable struct FFTData
  name::String                                # name [Net.Sta.Loc.Chan]
  id::String                                  # id [Y-mm-dd] this is date of fft
  loc::Array{Float64,1}                       # loc
  fs::Float64                                 # sampling rate [Hz]
  gain::Float64                               # gain
  freqmin::Float64                            # minumum frequency [Hz]
  freqmax::Float64                            # maximum frequency [Hz]
  cc_len::Int                                 # window_length [s]
  cc_step::Int                                # step between windows [s]
  whitened::Bool                              # whitening applied
  time_norm::Union{Bool,String}               # time normaliation
  resp  ::Array{Complex{Float64},2}           # response poles/zeros
  misc::Dict{String,Any}                      # misc
  notes::Array{String,1}                      # notes
  t::Array{Float64,1}                         # time
  fft::Array{Complex{Float64},2}              # fft data

  function FFTData(
      name     ::String,
      id       ::String,
      loc      ::Array{Float64,1},
      fs       ::Float64,
      gain     ::Float64,
      freqmin  ::Float64,
      freqmax  ::Float64,
      cc_len   ::Int,
      cc_step  ::Int,
      whitened ::Bool,
      time_norm::Union{Bool,String},
      resp     ::Array{Complex{Float64},2},
      misc     ::Dict{String,Any},
      notes    ::Array{String,1},
      t        ::Array{Float64,1},
      fft      ::Array{Complex{Float64},2}
      )

      return new(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
                 whitened, time_norm, resp, misc, notes, t, fft)
    end
end

FFTData(;
          name     ::String                    = "",
          id       ::String                    = "",
          loc      ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          fs       ::Float64                   = zero(Float64),
          gain     ::Float64                   = one(Float64),
          freqmin  ::Float64                   = zero(Float64),
          freqmax  ::Float64                   = zero(Float64),
          cc_len   ::Int                       = zero(Int),
          cc_step  ::Int                       = zero(Int),
          whitened ::Bool                      = false,
          time_norm::Union{Bool,String}        = false,
          resp     ::Array{Complex{Float64},2} = Array{Complex{Float64},2}(undef, 0, 2),
          misc     ::Dict{String,Any}          = Dict{String,Any}(),
          notes    ::Array{String,1}           = Array{String,1}(undef, 0),
          t        ::Array{Float64,1}            = Array{Float64,1}(undef, 0),
          fft      ::Array{Complex{Float64},2} = Array{Complex{Float64},2}(undef, 0, 2)
          ) = FFTData(name, id, loc, fs, gain, freqmin, freqmax, cc_len, cc_step,
                     whiten, time_norm, resp, misc, notes, t, fft)

FFTData(C::SeisChannel,freqmin::Float64, freqmax::Float64,cc_len::Int,
        cc_step::Int, whitened::Bool,time_norm::Union{Bool,String},
        t::Array{Float64,1},fft::Array{Complex{Float64},2}
       ) = FFTData(C.name, C.id, C.loc, C.fs, C.gain, freqmin, freqmax, cc_len,
                   cc_step, whitened, time_norm, C.resp, C.misc, C.notes, t, fft)

in(s::String, F::FFTData) = S.id==s

isempty(F::FFTData) = minimum([isempty(getfield(F,f)) for f in datafields])

isequal(F::FFTData, U::FFTData) = minimum([hash(getfield(S,i))==hash(getfield(U,i)) for i in fftfields]::Array{Bool,1})
==(F::FFTData, U::FFTData) = isequal(S,U)::Bool

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
