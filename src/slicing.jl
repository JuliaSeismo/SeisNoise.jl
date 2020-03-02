export start_end, slide, nearest_start_end, slide_ind
const μs = 1.0e-6
const sμ = 1000000.0
"""

  start_end(C::SeisChannel)

Return start and endtimes of SeisChannel in DateTime format.
"""
function start_end(C::SeisChannel)
    return u2d.(SeisIO.t_win(C.t,C.fs) * 1e-6)
end

"""

  start_end(S::SeisData)

Return start and endtimes of SeisChannel in DateTime format.
"""
function start_end(S::SeisData)
  irr = falses(S.n)
  non_ts = findall(S.fs .== 0)
  irr[non_ts] .= true
  z = zero(Int64)
  start_times = Array{DateTime}(undef,S.n)
  end_times =  Array{DateTime}(undef,S.n)
  for i = 1:S.n
    start_times[i] = u2d(S.t[i][1,2] * μs)
    end_times[i] = u2d(S.t[i][end,2] + (irr[i] ? z : S.t[i][1,2] + sum(S.t[i][2:end,2]) +
        round(Int64, (length(S.x[i])-1)/(μs*S.fs[i]))) * μs)
        # ts data: start time in μs from epoch + sum of time gaps in μs + length of trace in μs
        # non-ts data: time of last sample
  end
  return start_times,end_times
end

function slide(A::AbstractArray, cc_len::Int, cc_step::Int, fs::AbstractFloat,
               starttime::AbstractFloat,endtime::AbstractFloat)
    N = size(A,1)
    window_samples = Int(cc_len * fs)
    starts = Array(range(starttime,stop=endtime,step=cc_step))
    ends = starts .+ cc_len .- 1. / fs
    ind = findlast(x -> x <= endtime,ends)
    starts = starts[1:ind]
    ends = ends[1:ind]

    # fill array with overlapping windows
    if cc_step == cc_len && N % cc_len == 0
        return reshape(A,window_samples,N ÷ window_samples),starts
    elseif cc_step == cc_len # disregard data from edge
        return reshape(A[1 : N - N % window_samples], window_samples, N ÷ window_samples),starts
    else # need overlap between windows
        out = Array{eltype(A),2}(undef, window_samples,length(starts))
        s = convert.(Int,round.((hcat(starts,ends) .- starttime) .* fs .+ 1.))
        @inbounds for ii in eachindex(starts)
            out[:,ii] .= @view(A[s[ii,1]:s[ii,2]])
        end

    end
  return out,starts
end

"""
    slide(C::SeisChannel, cc_len::Int, cc_step::Int)

Cut `C` into equal length sliding windows.

# Arguments
- `C::SeisChannel`: SeisChannel.
- `cc_len::Int`: Cross-correlation window length [s].
- `cc_step::Int`: Step between cross-correlation windows [s].

# Returns
- `A::Array`: Array of sliding windows
- `starts::Array`: Array of start times of each window, in Unix time. E.g to convert
        Unix time to date time, use u2d(starts[1]) = 2018-08-12T00:00:00
- `ends::Array`: Array of end times of each window, in Unix time.
"""
slide(C::SeisChannel, cc_len::Int, cc_step::Int) = slide(C.x,cc_len,cc_step,C.fs,C.t)

"""
    nearest_start_end(C::SeisChannel, cc_len::Int, cc_step::Int)

Return best possible start, end times for data in `C` given the `cc_step` and `cc_len`.
"""
function nearest_start_end(C::SeisChannel, cc_len::Int, cc_step::Int)
  su,eu = SeisIO.t_win(C.t,C.fs) * μs
  ideal_start = d2u(DateTime(Date(u2d(su)))) # midnight of same day
  starts = Array(range(ideal_start,stop=eu,step=cc_step))
  ends = starts .+ cc_len .- 1. / C.fs
  return starts[findfirst(x -> x >= su, starts)], ends[findlast(x -> x <= eu,ends)]
end

"""
    nearest_start_end(S::DateTime,E::DateTime, cc_len::Int, cc_step::Int)

Return best possible start, end times for given starttime `S` and endtime `E`.
"""
function nearest_start_end(S::DateTime, E::DateTime, fs::Float64, cc_len::Int, cc_step::Int)
  ideal_start = DateTime(Date(S)) # midnight of same day
  starts = Array(ideal_start:Second(cc_step):endtime)
  ends = starts .+ Second(cc_len) .- Millisecond(convert(Int,1. / fs * 1e3))
  return starts[findfirst(x -> x >= S, starts)], ends[findlast(x -> x <= E,ends)]
end

function slide_ind(startslice::AbstractFloat,endslice::AbstractFloat,fs::AbstractFloat,t::AbstractArray)
  starttime = t[1,2] * 1e-6
  startind = convert(Int,round(startslice - starttime,digits=4) * fs) + 1
  endind = convert(Int,round(endslice - starttime,digits=4) * fs) + 1
  return startind,endind
end
