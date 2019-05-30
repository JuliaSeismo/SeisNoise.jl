import Base:show, size, summary
import Printf: @sprintf

const show_os = 12
const unindexed_fields = (:c, :n)
const FloatArray = Union{Array{Float64,1}, Array{Float32,1}}
# si(w::Int, i::Int) = show_os + w*(i-1)
# showtail(io::IO, b::Bool) = b ? "…" : ""
# float_str(x::Union{Float32,Float64}) = @sprintf("%.3e", x)
# # maxgap(t::Array{Int64,2}) = @sprintf("%g", μs*maximum(t[2:end,2]))
#
# function str_trunc(str::String, w::Integer)
#   d = map(UInt8, collect(str))
#   L = length(d)
#   if L > w
#     s3 = d[1:w-4]
#     d = map(UInt8, collect(String([s3; codeunits("...")])))
#   else
#     d = d[1:min(w,L)]
#   end
#   return d
# end
#
# function str_head(s::String, W::Int)
#   sd = ones(UInt8, W)*0x20
#   sd[show_os-1-length(s):show_os-2] = map(UInt8, collect(uppercase(s)))
#   sd[show_os-1:show_os] = codeunits(": ")
#   return sd
# end
#
# function show_str(io::IO, S::Array{String,1}, w::Int, W::Int, s::String, b::Bool)
#   sd = str_head(s, W)
#   N = length(S)
#   for i = 1:N
#     st = si(w,i)
#     d = str_trunc(S[i], w)
#     sd[st+1:st+length(d)] = d
#   end
#   println(io, replace(String(sd),'\0' => ' '), showtail(io, b))
#   return
# end
#
# function show_t(io::IO, T::Array{Float64,1}, w::Int, W::Int, b::Bool)
#   sd1 = str_head("T", W::Int)
#   p = show_os
#   for i = 1:length(T)
#     if isempty(T[i])
#       s = ""
#     else
#       s = u2d.(T[i])
#     end
#     sd1[p+1:p+length(s)] = codeunits(s)
#     p += w
#   end
#   println(io, replace(String(sd1),'\0' => ' '), showtail(io, b))
#   return
# end
#
# function show_x(io::IO,
#                 X::Union{ Array{Array{Float64,1},1},
#                           Array{Array{Float32,1},1},
#                           Array{Union{Array{Float64,1}, Array{Float32,1}},1} },
#                 w::Int, W::Int, tip::String, b::Bool)
#   N = length(X)
#   str = zeros(UInt8, W, 6)
#   str[show_os-length(tip)-1:show_os,1] = UInt8.(codeunits(tip * ": "))
#   p = show_os
#   i = 1
#   while p < W && i <= N
#     L = length(X[i])
#     Lx = string("(nx = ", L, ")")
#     if isempty(X[i])
#       str[p+1:p+7,1] = codeunits("(empty)")
#     else
#       for k = 1:5
#         if k <= L
#           s = float_str(X[i][k])
#           if (L > 5 && k==3)
#             s = "  ..."
#           elseif (L > 5 && k==4)
#             s = float_str(last(X[i]))
#           elseif (L > 5 && k==5)
#             s = Lx
#           end
#         else
#           s = ""
#         end
#         cstr = codeunits(s)
#         str[p+1:p+length(cstr),k] = cstr
#       end
#     end
#     p += w
#     i += 1
#   end
#   for i = 1:5
#     if i == 1
#       println(io, replace(String(str[:,i]),'\0' => ' '), showtail(io, b))
#     else
#       println(io, replace(String(str[:,i]),'\0' => ' '))
#     end
#   end
#   return
# end

# function resp_str(io::IO, X::Array{Array{Complex{Float64},2},1}, w::Int, W::Int, b::Bool)
#   N = length(X); p = show_os; i = 1
#   sd = zeros(UInt8, W, 2)
#   sd[show_os-5:show_os-1,1] = codeunits("RESP:")
#   while p < W && i <= N
#     zstr = ""
#     pstr = ""
#     R = X[i]
#     if isempty(R)
#       zstr *= " "
#       pstr *= " "
#     else
#       L = size(R,1)
#       for j = 1:L
#         zstr *= string(float_str(real(R[j,1])), " + ", float_str(imag(R[j,1])), "i")
#         pstr *= string(float_str(real(R[j,2])), " + ", float_str(imag(R[j,2])), "i")
#         if L > 1 && j < L
#           zstr *= ", "
#           pstr *= ", "
#         end
#         if length(zstr) > w
#           zstr = zstr[1:w-4]*"..."
#           pstr = pstr[1:w-4]*"..."
#           break
#         end
#       end
#     end
#     q = length(zstr)
#     sd[p+1:p+q,1] = codeunits(zstr)
#     sd[p+1:p+q,2] = codeunits(pstr)
#     p += w
#     i += 1
#   end
#   [println(io, replace(String(sd[:,i]),'\0' => ' '), showtail(io,b)) for i=1:2]
#   return
# end

showtail(b::Bool) = b ? "…" : ""
ngaps(t::Array{Int64,2}) = max(0, size(t,1)-2)

function str_trunc(s::String, w::Int64)
  L = length(s)
  if L ≥ w
    return s[1:w-2] * "… "
  else
    return rpad(s, w)
  end
end

function show_str(io::IO, S::Array{String,1}, w::Int, N::Int64)
  for i = 1:N
    print(io, str_trunc(S[i], w))
  end
  println(io, showtail(N<length(S)))
  return nothing
end

function show_float(io::IO, X::FloatArray, w::Int, N::Int64)
  for i = 1:N
     print(io, rpad(@sprintf("%+10.3e", X[i]), w))
  end
  println(io, showtail(N<length(S)))
  return nothing
end

function show_t(io::IO, T::Array{Array{Int64,2},1}, w::Int, N::Int64)
  for i = 1:N
    if isempty(T[i])
      s = ""
    else
      s = string(timestamp(T[i][1,2]*μs), " (", ngaps(T[i]), " gaps)")
    end
    print(io, str_trunc(s, w))
  end
  println(io, showtail(N<length(T)))
  return
end

function mkxstr(N, X::Union{Array{Array{Float64,1},1},
                            Array{Array{Float32,1},1},
                            Array{Union{Array{Float64,1}, Array{Float32,1}},1}})

  # Fill matrix of X values
  vx = 5
  X_str = Array{String,2}(undef, vx, N)
  fill!(X_str, "")
  for j = 1:N
    x = getindex(X, j)
    nx = lastindex(x)
    if nx == 0
      X_str[1,j] = "(empty)"
      continue
    elseif nx < vx
      for i = 1:nx
        X_str[i,j] = @sprintf("%+10.3e", x[i])
      end
    else
      nx_str          = string(nx)
      for i = 1:vx-3
        X_str[i,j]    = @sprintf("%+10.3e", x[i])
      end
      X_str[vx-2,j]   = "    ..."
      X_str[vx-1,j]   = @sprintf("%+10.3e", x[nx])
      X_str[vx,j]     = string("(nx = ", nx_str, ")")
    end
  end
  return X_str
end

# Fill matrix of X value strings
function mkxstr(X::FloatArray)
  vx = 5
  X_str = Array{String,2}(undef, vx, 1)
  fill!(X_str, "")
  nx = lastindex(X)
  if nx == 0
    X_str[1,1]      = "(empty)"
  elseif nx < vx
    for i = 1:nx
      X_str[i,1]    = @sprintf("%+10.3e", X[i])
    end
  else
    nx_str          = string(nx)
    for i = 1:vx-3
      X_str[i,1]    = @sprintf("%+10.3e", X[i])
    end
    X_str[vx-2,1]   = "    ..."
    X_str[vx-1,1]   = @sprintf("%+10.3e", X[nx])
    X_str[vx,1]     = string("(nx = ", nx_str, ")")
  end
  return X_str
end

function show_x(io::IO, X_str::Array{String,2}, w::Int64, b::Bool)
  (vx, N) = size(X_str)

  # Display
  for i = 1:vx
    if i > 1
      print(io, " "^show_os)
    end

    for j = 1:N
      x_str = X_str[i,j]
      L = length(x_str)
      print(io, x_str)
      if (x_str == "(empty)" || x_str == "") && N == 1
        return nothing
      end
      print(io, " "^(w-L))
    end
    print(io, showtail(b))
    if i < vx
      print(io, "\n")
    end
  end
  return nothing
end

function showloc_full(io::IO, Loc::T) where {T<:InstrumentPosition}
  F = fieldnames(T)
  println(io, T, " with fields:")
  for f in F
    fn = lpad(String(f), 5, " ")
    println(io, fn, ": ", getfield(Loc,f))
  end
  return nothing
end

function show(io::IO, Loc::GenLoc)
  if get(io, :compact, false) == false
    showloc_full(io, Loc)
  else
    print(io, repr(getfield(Loc, :loc), context=:compact => true))
  end
  return nothing
end

function show(io::IO, loc::GeoLoc)
  if get(io, :compact, false) == false
    showloc_full(io, loc)
  else
    c = :compact => true
    print(io, repr(getfield(loc, :lat), context=c), " N, ",
              repr(getfield(loc, :lon), context=c), " E, ",
              repr(getfield(loc, :el), context=c), " m")
  end
  return nothing
end


function show(io::IO, loc::UTMLoc)
  if get(io, :compact, false) == false
    showloc_full(io, loc)
  else
    print(io, getfield(loc, :zone), " ",
              getfield(loc, :hemi), " ",
              getfield(loc, :E), " ",
              getfield(loc, :N))
  end
  return nothing
end

function show(io::IO, loc::XYLoc)
  if get(io, :compact, false) == false
    showloc_full(io, loc)
  else
    c = :compact => true
    print(io, "x ", repr(getfield(loc, :x), context=c),
              ", y ", repr(getfield(loc, :x), context=c),
              ", z ", repr(getfield(loc, :x), context=c))
  end
  return nothing
end

function showresp_full(io::IO, Resp::T) where {T<:InstrumentResponse}
  F = fieldnames(T)
  println(io, T, " with fields:")
  for f in F
    fn = lpad(String(f), 5, " ")
    println(io, fn, ": ", getfield(Resp,f))
  end
  return nothing
end

function show(io::IO, Resp::GenResp)
  if get(io, :compact, false) == false
    showresp_full(io, Resp)
  else
    resp = getfield(Resp, :resp)
    M,N = size(resp)
    M1 = min(M,2)
    N1 = min(N,2)
    print(io, "[")
    for i = 1:M1
      for j = 1:N1
        print(io, repr(resp[i,j], context=:compact=>true))
        if j == N1 && i < M1
          if N > N1
            print(io, " … ")
          end
          print(io, "; ")
        elseif j == N1 && i == M1
          if N > N1
            print(io, " … ;")
          end
          if M > M1
            print(io, " … ")
          end
          print(io, "]")
        else
          print(io, ", ")
        end
      end
    end
    print(io, " (")
    print(io, getfield(Resp, :desc))
    print(io, ")")
  end
  return nothing
end

function show(io::IO, Resp::Union{PZResp,PZResp64})
  if get(io, :compact, false) == false
    showresp_full(io, Resp)
  else
    c = :compact => true
    print(io, "c = ", repr(getfield(Resp, :c), context=c), ", ",
              length(getfield(Resp, :z)), " zeros, ",
              length(getfield(Resp, :p)), " poles")
  end
  return nothing
end

summary(F::FFTData) = string(typeof(F), " with ", length(F.t), " fft",
                             size(F.fft,2) == 1 ? "" : "s")
summary(C::CorrData) = string(typeof(C), " with ", length(C.t),
                              " correlation", size(C.corr,1) == 1 ? "" : "s")

# function show(io::IO, F::FFTData)
#   loc_str = ["lat", "lon", "ele", "az", "inc"]
#   W = max(80,displaysize(io)[2]-2)-show_os
#   w = min(W, 36)
#   D = Array{String,1}(undef,25)
#   fftlen,numfft = size(F.fft)
#   println(io, "FFTData with ", numfft, " ffts")
#   show_str(io, [F.id], w, W, "id", false)
#   show_str(io, [F.name], w, W, "name", false)
#   show_x(io, [F.loc], w, W, "LOC", false)
#   show_str(io, [string(F.fs)], w, W, "fs", false)
#   show_str(io, [string(F.gain)], w, W, "gain", false)
#   show_str(io, [string(F.freqmin)], w, W, "freqmin", false)
#   show_str(io, [string(F.freqmax)], w, W, "freqmax", false)
#   show_str(io, [string(F.cc_len)], w, W, "cc_len", false)
#   show_str(io, [string(F.cc_step)], w, W, "cc_step", false)
#   show_str(io, [string(F.whitened)], w, W, "whitened", false)
#   show_str(io, [string(F.time_norm)], w, W, "cc_norm", false)
#   resp_str(io, [F.resp], w, W, false)
#   show_str(io, [string(length(F.notes), " entries")], w, W, "NOTES", false)
#   show_str(io, [string(length(F.misc), " items")], w, W, "MISC", false)
#   # show_t(io, [F.t], w, W, false)
#   show_str(io, [string(u2d(F.t[1]))], w, W, "starttime",false)
#   show_str(io, [string((u2d(F.t[end] + F.cc_len - 1 / F.fs)))],
#                 w, W, "endtime",false)
#   # show(io, "text/plain", F.fft)
#   return nothing
# end
# show(F::FFTData) = show(stdout, F)
#
# function show(io::IO, C::CorrData)
#   loc_str = ["lat", "lon", "ele", "az", "inc"]
#   W = max(80,displaysize(io)[2]-2)-show_os
#   w = min(W, 36)
#   D = Array{String,1}(undef,25)
#   corrlen,numcorr = size(C.corr)
#   println(io, "CorrData with ", numcorr, " ffts")
#   show_str(io, [C.id], w, W, "id", false)
#   show_str(io, [C.name], w, W, "name", false)
#   show_x(io, [C.loc], w, W, "LOC", false)
#   show_str(io, [C.comp], w, W, "comp", false)
#   show_str(io, [string(C.rotated)], w, W, "rotated", false)
#   show_str(io, [string(C.corr_type)], w, W, "corr_type", false)
#   show_str(io, [string(C.fs)], w, W, "fs", false)
#   show_str(io, [string(C.gain)], w, W, "gain", false)
#   show_str(io, [string(C.freqmin)], w, W, "freqmin", false)
#   show_str(io, [string(C.freqmax)], w, W, "freqmax", false)
#   show_str(io, [string(C.cc_len)], w, W, "cc_len", false)
#   show_str(io, [string(C.cc_step)], w, W, "cc_step", false)
#   show_str(io, [string(C.whitened)], w, W, "whitened", false)
#   show_str(io, [string(C.time_norm)], w, W, "cc_norm", false)
#   resp_str(io, [C.resp], w, W, false)
#   show_str(io, [string(length(C.notes), " entries")], w, W, "NOTES", false)
#   show_str(io, [string(length(C.misc), " items")], w, W, "MISC", false)
#   show_str(io, [string(C.maxlag)], w, W, "maxlag", false)
#   show_str(io, [string(u2d(C.t[1]))], w, W, "starttime",false)
#   show_str(io, [string((u2d(C.t[end] + C.cc_len - 1 / C.fs)))],
#                 w, W, "endtime",false)
#   return nothing
# end
# show(C::CorrData) = show(stdout, C)

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
      if t == Array{String,1}
        SeisIO.show_str(io, targ, w, N)
      elseif f == :notes || f == :misc
        SeisIO.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(F.t) > 0
          SeisIO.show_str(io, String[timestamp(F.t[1]), " (", string(length(F.t)), " FFTs)"], w, N)
        else
          SeisIO.show_str(io, String[repr("text/plain", targ[i], context=:compact=>true) for i = 1:M], w, N)
        end
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
      if t == Array{String,1}
        SeisIO.show_str(io, targ, w, N)
      elseif f == :notes || f == :misc
        SeisIO.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(C.t) > 0
          SeisIO.show_str(io, String[timestamp(C.t[1]), " (", string(length(C.t)), " Corrs)"], w, N)
        else
          SeisIO.show_str(io, String[repr("text/plain", targ[i], context=:compact=>true) for i = 1:M], w, N)
        end
      else
        SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
      end
    end
  end
  return nothing
end
show(C::CorrData) = show(stdout, C)
