import Base:show, size, summary
import Printf: @sprintf

const show_os = 12
si(w::Int, i::Int) = show_os + w*(i-1)
showtail(io::IO, b::Bool) = b ? "…" : ""
float_str(x::Union{Float32,Float64}) = @sprintf("%.3e", x)
# maxgap(t::Array{Int64,2}) = @sprintf("%g", μs*maximum(t[2:end,2]))

function str_trunc(str::String, w::Integer)
  d = map(UInt8, collect(str))
  L = length(d)
  if L > w
    s3 = d[1:w-4]
    d = map(UInt8, collect(String([s3; codeunits("...")])))
  else
    d = d[1:min(w,L)]
  end
  return d
end

function str_head(s::String, W::Int)
  sd = ones(UInt8, W)*0x20
  sd[show_os-1-length(s):show_os-2] = map(UInt8, collect(uppercase(s)))
  sd[show_os-1:show_os] = codeunits(": ")
  return sd
end

function show_str(io::IO, S::Array{String,1}, w::Int, W::Int, s::String, b::Bool)
  sd = str_head(s, W)
  N = length(S)
  for i = 1:N
    st = si(w,i)
    d = str_trunc(S[i], w)
    sd[st+1:st+length(d)] = d
  end
  println(io, replace(String(sd),'\0' => ' '), showtail(io, b))
  return
end

function show_t(io::IO, T::Array{Float64,1}, w::Int, W::Int, b::Bool)
  sd1 = str_head("T", W::Int)
  p = show_os
  for i = 1:length(T)
    if isempty(T[i])
      s = ""
    else
      s = u2d.(T[i])
    end
    sd1[p+1:p+length(s)] = codeunits(s)
    p += w
  end
  println(io, replace(String(sd1),'\0' => ' '), showtail(io, b))
  return
end

function show_x(io::IO,
                X::Union{ Array{Array{Float64,1},1},
                          Array{Array{Float32,1},1},
                          Array{Union{Array{Float64,1}, Array{Float32,1}},1} },
                w::Int, W::Int, tip::String, b::Bool)
  N = length(X)
  str = zeros(UInt8, W, 6)
  str[show_os-length(tip)-1:show_os,1] = UInt8.(codeunits(tip * ": "))
  p = show_os
  i = 1
  while p < W && i <= N
    L = length(X[i])
    Lx = string("(nx = ", L, ")")
    if isempty(X[i])
      str[p+1:p+7,1] = codeunits("(empty)")
    else
      for k = 1:5
        if k <= L
          s = float_str(X[i][k])
          if (L > 5 && k==3)
            s = "  ..."
          elseif (L > 5 && k==4)
            s = float_str(last(X[i]))
          elseif (L > 5 && k==5)
            s = Lx
          end
        else
          s = ""
        end
        cstr = codeunits(s)
        str[p+1:p+length(cstr),k] = cstr
      end
    end
    p += w
    i += 1
  end
  for i = 1:5
    if i == 1
      println(io, replace(String(str[:,i]),'\0' => ' '), showtail(io, b))
    else
      println(io, replace(String(str[:,i]),'\0' => ' '))
    end
  end
  return
end

function resp_str(io::IO, X::Array{Array{Complex{Float64},2},1}, w::Int, W::Int, b::Bool)
  N = length(X); p = show_os; i = 1
  sd = zeros(UInt8, W, 2)
  sd[show_os-5:show_os-1,1] = codeunits("RESP:")
  while p < W && i <= N
    zstr = ""
    pstr = ""
    R = X[i]
    if isempty(R)
      zstr *= " "
      pstr *= " "
    else
      L = size(R,1)
      for j = 1:L
        zstr *= string(float_str(real(R[j,1])), " + ", float_str(imag(R[j,1])), "i")
        pstr *= string(float_str(real(R[j,2])), " + ", float_str(imag(R[j,2])), "i")
        if L > 1 && j < L
          zstr *= ", "
          pstr *= ", "
        end
        if length(zstr) > w
          zstr = zstr[1:w-4]*"..."
          pstr = pstr[1:w-4]*"..."
          break
        end
      end
    end
    q = length(zstr)
    sd[p+1:p+q,1] = codeunits(zstr)
    sd[p+1:p+q,2] = codeunits(pstr)
    p += w
    i += 1
  end
  [println(io, replace(String(sd[:,i]),'\0' => ' '), showtail(io,b)) for i=1:2]
  return
end


summary(F::FFTData) = string(typeof(F), " with ", size(F.fft,2), " fft",
                             size(F.fft,2) == 1 ? "" : "s")
summary(C::CorrData) = string(typeof(C), " with ", size(C.corr,2),
                              " correlation", size(C.corr,2) == 1 ? "" : "s")

function show(io::IO, F::FFTData)
  loc_str = ["lat", "lon", "ele", "az", "inc"]
  W = max(80,displaysize(io)[2]-2)-show_os
  w = min(W, 36)
  D = Array{String,1}(undef,25)
  fftlen,numfft = size(F.fft)
  println(io, "FFTData with ", numfft, " ffts")
  show_str(io, [F.id], w, W, "id", false)
  show_str(io, [F.name], w, W, "name", false)
  show_x(io, [F.loc], w, W, "LOC", false)
  show_str(io, [string(F.fs)], w, W, "fs", false)
  show_str(io, [string(F.gain)], w, W, "gain", false)
  show_str(io, [string(F.freqmin)], w, W, "freqmin", false)
  show_str(io, [string(F.freqmax)], w, W, "freqmax", false)
  show_str(io, [string(F.cc_len)], w, W, "cc_len", false)
  show_str(io, [string(F.cc_step)], w, W, "cc_step", false)
  show_str(io, [string(F.whitened)], w, W, "whitened", false)
  show_str(io, [string(F.time_norm)], w, W, "cc_norm", false)
  resp_str(io, [F.resp], w, W, false)
  show_str(io, [string(length(F.notes), " entries")], w, W, "NOTES", false)
  show_str(io, [string(length(F.misc), " items")], w, W, "MISC", false)
  # show_t(io, [F.t], w, W, false)
  show_str(io, [string(u2d(F.t[1]))], w, W, "starttime",false)
  show_str(io, [string((u2d(F.t[end] + F.cc_len - 1 / F.fs)))],
                w, W, "endtime",false)
  # show(io, "text/plain", F.fft)
  return nothing
end
show(F::FFTData) = show(stdout, F)

function show(io::IO, C::CorrData)
  loc_str = ["lat", "lon", "ele", "az", "inc"]
  W = max(80,displaysize(io)[2]-2)-show_os
  w = min(W, 36)
  D = Array{String,1}(undef,25)
  corrlen,numcorr = size(C.corr)
  println(io, "CorrData with ", numcorr, " ffts")
  show_str(io, [C.id], w, W, "id", false)
  show_str(io, [C.name], w, W, "name", false)
  show_x(io, [C.loc], w, W, "LOC", false)
  show_str(io, [C.comp], w, W, "comp", false)
  show_str(io, [string(C.rotated)], w, W, "rotated", false)
  show_str(io, [string(C.corr_type)], w, W, "corr_type", false)
  show_str(io, [string(C.fs)], w, W, "fs", false)
  show_str(io, [string(C.gain)], w, W, "gain", false)
  show_str(io, [string(C.freqmin)], w, W, "freqmin", false)
  show_str(io, [string(C.freqmax)], w, W, "freqmax", false)
  show_str(io, [string(C.cc_len)], w, W, "cc_len", false)
  show_str(io, [string(C.cc_step)], w, W, "cc_step", false)
  show_str(io, [string(C.whitened)], w, W, "whitened", false)
  show_str(io, [string(C.time_norm)], w, W, "cc_norm", false)
  resp_str(io, [C.resp], w, W, false)
  show_str(io, [string(length(C.notes), " entries")], w, W, "NOTES", false)
  show_str(io, [string(length(C.misc), " items")], w, W, "MISC", false)
  show_str(io, [string(C.maxlag)], w, W, "maxlag", false)
  show_str(io, [string(u2d(C.t[1]))], w, W, "starttime",false)
  show_str(io, [string((u2d(C.t[end] + C.cc_len - 1 / C.fs)))],
                w, W, "endtime",false)
  return nothing
end
show(C::CorrData) = show(stdout, C)
