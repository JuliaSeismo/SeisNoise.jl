import Base:show, size, summary

const show_os = 12
const unindexed_fields = (:c, :n)

summary(R::RawData) = string(typeof(R), " with ", length(R.t),
                              " window", size(R.x,1) == 1 ? "" : "s")
summary(F::FFTData) = string(typeof(F), " with ", length(F.t), " fft",
                             size(F.fft,2) == 1 ? "" : "s")
summary(C::CorrData) = string(typeof(C), " with ", length(C.t),
                              " correlation", size(C.corr,1) == 1 ? "" : "s")

function show(io::IO, R::RawData)
  W = max(80, displaysize(io)[2]) - show_os
  w = min(W, 35)
  nc = size(R.x,1) == 0 ? 0 : 1
  N = min(nc, div(W-1, w))
  M = min(N+1, nc)
  rawlen,numraw = size(R.x)
  println(io, "RawData with ", length(R.t), " windows")
  Field = fieldnames(RawData)
  for f in Field
    if (f in unindexed_fields) == false
      targ = getfield(R, f)
      t = typeof(targ)
      fstr = uppercase(String(f))
      print(io, lpad(fstr, show_os-2), ": ")
      if f == :notes || f == :misc
        SeisIO.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(R.t) > 0
          SeisIO.show_str(io, String[timestamp(R.t[1]), " (", string(length(R.t)), " windows)"], w, N)
        else
          SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
        end
      elseif f == :x
        SeisIO.show_str(io, String[summary(targ)], w, N)
      else
        SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
      end
    end
  end
  return nothing
end
show(R::RawData) = show(stdout, R)


function show(io::IO, F::FFTData)
  W = max(80, displaysize(io)[2]) - show_os
  w = min(W, 35)
  nc = size(F.fft,1) == 0 ? 0 : 1
  N = min(nc, div(W-1, w))
  M = min(N+1, nc)
  fftlen,numfft = size(F.fft)
  println(io, "FFTData with ", length(F.t), " ffts")
  Field = fieldnames(FFTData)
  for f in Field
    if (f in unindexed_fields) == false
      targ = getfield(F, f)
      t = typeof(targ)
      fstr = uppercase(String(f))
      print(io, lpad(fstr, show_os-2), ": ")
      if f == :notes || f == :misc
        SeisIO.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(F.t) > 0
          SeisIO.show_str(io, String[timestamp(F.t[1]), " (", string(length(F.t)), " FFTs)"], w, N)
        else
          SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
        end
      elseif f == :fft
        SeisIO.show_str(io, String[summary(targ)], w, N)
      else
        SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
      end
    end
  end
  return nothing
end
show(F::FFTData) = show(stdout, F)

function show(io::IO, C::CorrData)
  W = max(80, displaysize(io)[2]) - show_os
  w = min(W, 35)
  nc = size(C.corr,1) == 0 ? 0 : 1
  N = min(nc, div(W-1, w))
  M = min(N+1, nc)
  corrlen,numcorr = size(C.corr)
  println(io, "CorrData with ", length(C.t), " Corrs")
  Field = fieldnames(CorrData)
  for f in Field
    if (f in unindexed_fields) == false
      targ = getfield(C, f)
      t = typeof(targ)
      fstr = uppercase(String(f))
      print(io, lpad(fstr, show_os-2), ": ")
      if f == :notes || f == :misc
        SeisIO.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(C.t) > 0
          SeisIO.show_str(io, String[timestamp(C.t[1]), " (", string(length(C.t)), " Corrs)"], w, N)
        else
          SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
        end
      elseif f == :corr
        SeisIO.show_str(io, String[summary(targ)], w, N)
      else
        SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
      end
    end
  end
  return nothing
end
show(C::CorrData) = show(stdout, C)
