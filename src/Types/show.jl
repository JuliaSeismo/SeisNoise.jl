import Base:show, size, summary

const show_os = 12
const unindexed_fields = (:n, :ns, :ox, :oy, :oz, :freqmax ,:freqmin, :cc_len, :maxlag, :whitened,
                          :time_norm, :data, :info, :x, :preprocessed, :dims)

summary(R::RawData) = string(typeof(R), " with ", length(R.t),
                              " window", size(R.x,1) == 1 ? "" : "s")
summary(F::FFTData) = string(typeof(F), " with ", length(F.t), " fft",
                             size(F.fft,2) == 1 ? "" : "s")
summary(C::CorrData) = string(typeof(C), " with ", length(C.t),
                              " correlation", size(C.corr,1) == 1 ? "" : "s")
# need to add summary functions for NodalFFTData, NodalProcessedData, and NodalCorrData

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
        SeisBase.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(R.t) > 0
          SeisBase.show_str(io, String[timestamp(R.t[1]), " (", string(length(R.t)), " windows)"], w, N)
        else
          SeisBase.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
        end
      elseif f == :x
        SeisBase.show_str(io, String[summary(targ)], w, N)
      else
        SeisBase.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
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
        SeisBase.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(F.t) > 0
          SeisBase.show_str(io, String[timestamp(F.t[1]), " (", string(length(F.t)), " FFTs)"], w, N)
        else
          SeisBase.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
        end
      elseif f == :fft
        SeisBase.show_str(io, String[summary(targ)], w, N)
      else
        SeisBase.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
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
        SeisBase.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(C.t) > 0
          SeisBase.show_str(io, String[timestamp(C.t[1]), " (", string(length(C.t)), " Corrs)"], w, N)
        else
          SeisBase.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
        end
      elseif f == :corr
        SeisBase.show_str(io, String[summary(targ)], w, N)
      else
        SeisBase.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
      end
    end
  end
  return nothing
end
show(C::CorrData) = show(stdout, C)

function show(io::IO, NF::NodalFFTData)
    W = max(80, displaysize(io)[2]) - show_os
    nc = getfield(NF, :n)
    w = min(W, 35)
    N = min(nc, div(W-1, w))
    M = min(N+1, nc)
    println(io, "NodalFFTData with ", nc, " channels (", N, " shown)")
    F = fieldnames(NodalFFTData)
    for f in F
        if ((f in unindexed_fields) == false) || (f == :x)
            targ = getfield(NF, f)
            t = typeof(targ)
            fstr = uppercase(String(f))
            print(io, lpad(fstr, show_os-2), ": ")
            if t == Array{String,1}
                show_str(io, targ, w, N)
            elseif f == :notes || f == :misc
                show_str(io, String[string(length(getindex(targ, i)), " entries") for i = 1:M], w, N)
            elseif f == :t
                show_t(io, targ, w, N, NF.fs)
            elseif f == :x
                x_str = mkxstr(N, getfield(NF, :x))
                show_x(io, x_str, w, N<nc)
            else
                show_str(io, String[repr("text/plain", targ[i], context=:compact=>true) for i = 1:M], w, N)
            end
        elseif f == :ox
            print(io, "COORDS: X = ", repr("text/plain", getfield(NF, f), context=:compact=>true), ", ")
        elseif f == :oy
            print(io, "Y = ", repr("text/plain", getfield(NF, f), context=:compact=>true), ", ")
        elseif f == :oz
            print(io, "Z = ", repr("text/plain", getfield(NF, f), context=:compact=>true), "\n")
        elseif f == :info
            print(io, "  INFO: ", length(NF.info), " entries\n")
        elseif f == :dims
            if NF.dims == [1]
                print(io, "  DIMS: time\n")
            elseif NF.dims == [2]
                print(io, "  DIMS: space\n")
            elseif NF.dims == [1,2]
                print(io, "  DIMS: time, space\n")
            end
        end
    end
    return nothing
end
show(NF::NodalFFTData) = show(stdout, NF)

function show(io::IO, NP::NodalProcessedData)
    W = max(80, displaysize(io)[2]) - show_os
    nc = getfield(NP, :n)
    w = min(W, 35)
    N = min(nc, div(W-1, w))
    M = min(N+1, nc)
    print()
    println(io, "NodalProcessedData with ", nc, " channels (", N, " shown)")
    F = fieldnames(NodalProcessedData)
    for f in F
        if ((f in unindexed_fields) == false) || (f == :x)
            targ = getfield(NP, f)
            t = typeof(targ)
            fstr = uppercase(String(f))
            print(io, lpad(fstr, show_os-2), ": ")
            if t == Array{String,1}
                show_str(io, targ, w, N)
            elseif f == :notes || f == :misc
                show_str(io, String[string(length(getindex(targ, i)), " entries") for i = 1:M], w, N)
            elseif f == :t
                show_t(io, targ, w, N, NP.fs)
            elseif f == :x
                x_str = mkxstr(N, getfield(NP, :x))
                show_x(io, x_str, w, N<nc)
            else
                show_str(io, String[repr("text/plain", targ[i], context=:compact=>true) for i = 1:M], w, N)
            end
        elseif f == :ox
            print(io, "COORDS: X = ", repr("text/plain", getfield(NP, f), context=:compact=>true), ", ")
        elseif f == :oy
            print(io, "Y = ", repr("text/plain", getfield(NP, f), context=:compact=>true), ", ")
        elseif f == :oz
            print(io, "Z = ", repr("text/plain", getfield(NP, f), context=:compact=>true), "\n")
        elseif f == :info
            print(io, "  INFO: ", length(NP.info), " entries\n")
        end
    end
    return nothing
end
show(NP::NodalProcessedData) = show(stdout, NP)

function show(io::IO, NC::NodalCorrData)
    W = max(80, displaysize(io)[2]) - show_os
    nc = getfield(NC, :n)
    w = min(W, 35)
    N = min(nc, div(W-1, w))
    M = min(N+1, nc)
    println(io, "NodalCorrData with ", nc, " correlations (", N, " shown)")
    F = fieldnames(NodalCorrData)
    for f in F
        if ((f in unindexed_fields) == false) || (f == :x)
            targ = getfield(NC, f)
            t = typeof(targ)
            fstr = uppercase(String(f))
            print(io, lpad(fstr, show_os-2), ": ")
            if t == Array{String,1}
                show_str(io, targ, w, N)
            elseif f == :notes || f == :misc
                show_str(io, String[string(length(getindex(targ, i)), " entries") for i = 1:M], w, N)
            elseif f == :t
                show_t(io, targ, w, N, NC.fs)
            elseif f == :x
                x_str = mkxstr(N, getfield(NC, :x))
                show_x(io, x_str, w, N<nc)
            else
                show_str(io, String[repr("text/plain", targ[i], context=:compact=>true) for i = 1:M], w, N)
            end
        elseif f == :ox
            print(io, "COORDS: X = ", repr("text/plain", getfield(NC, f), context=:compact=>true), ", ")
        elseif f == :oy
            print(io, "Y = ", repr("text/plain", getfield(NC, f), context=:compact=>true), ", ")
        elseif f == :oz
            print(io, "Z = ", repr("text/plain", getfield(NC, f), context=:compact=>true), "\n")
        elseif f == :info
            print(io, "  INFO: ", length(NC.info), " entries\n")
        end
    end
    return nothing
end
show(NC::NodalCorrData) = show(stdout, NC)
