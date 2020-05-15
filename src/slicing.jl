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

"""
    slide(A, cc_len, cc_step, fs, starttime, endtime)

Cut `A` into sliding windows of length `cc_len` points and offset `cc_step` points.

# Arguments
- `A::AbstractArray`: 1D time series.
- `cc_len::Float64`: Cross-correlation window length in seconds.
- `cc_step::Float64`: Step between cross-correlation windows in seconds.
- `starttime::Float64`: Time of first sample in `A` in unix time.
- `endtime::Float64`: Time of last sample in `A` in unix time.

# Returns
- `out::Array`: Array of sliding windows
- `starts::Array`: Array of start times of each window, in Unix time. E.g to convert
        Unix time to date time, use u2d(starts[1]) = 2018-08-12T00:00:00
"""
function slide(A::AbstractArray, cc_len::Float64, cc_step::Float64, fs::AbstractFloat,
               starttime::Float64,endtime::Float64)
    N = size(A,1)
    window_samples = Int(cc_len * fs)
    starts = Array(range(starttime,stop=endtime,step=cc_step))
    ends = starts .+ cc_len .- 1. / fs
    ind = findlast(x -> x <= endtime,ends)
    starts = starts[1:ind]
    ends = ends[1:ind]

    # fill array with overlapping windows
    if cc_step == cc_len && N % cc_len == 0
        return Array(reshape(A,window_samples,N ÷ window_samples)),starts
    elseif cc_step == cc_len # disregard data from edge
        return Array(reshape(A[1 : N - N % window_samples], window_samples, N ÷ window_samples)),starts
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
    slide(C, cc_len, cc_step)

Cut `C` into equal length sliding windows.

# Arguments
- `C::SeisChannel`: SeisChannel.
- `cc_len::Float64`: Cross-correlation window length in seconds.
- `cc_step::Float64`: Step between cross-correlation windows in seconds.

# Returns
- `out::Array`: Array of sliding windows.
- `starts::Array`: Array of start times of each window, in Unix time. E.g to convert
        Unix time to date time, use u2d(starts[1]) = 2018-08-12T00:00:00
"""
slide(C::SeisChannel, cc_len::Float64, cc_step::Float64) = slide(C.x,cc_len,cc_step,C.fs,d2u.(start_end(C))...)

"""
    nearest_start_end(C::SeisChannel, cc_len::Float64, cc_step::Float64)

Return best possible start, end times for data in `C` given the `cc_step` and `cc_len`.
"""
function nearest_start_end(C::SeisChannel, cc_len::Float64, cc_step::Float64)
  su,eu = SeisIO.t_win(C.t,C.fs) * μs
  eu = round(eu,digits=4) # round due to numerical roundoff error
  ideal_start = d2u(DateTime(Date(u2d(su)))) # midnight of same day
  starts = Array(range(ideal_start,stop=eu,step=cc_step))
  ends = starts .+ cc_len .- 1. / C.fs
  return starts[findfirst(x -> x >= su, starts)], ends[findlast(x -> x <= eu,ends)]
end

"""
    nearest_start_end(S::DateTime,E::DateTime, cc_len::Int, cc_step::Int)

Return best possible start, end times for given starttime `S` and endtime `E`.
"""
function nearest_start_end(S::DateTime, E::DateTime, fs::Float64, cc_len::Float64, cc_step::Float64)
  ideal_start = DateTime(Date(S)) # midnight of same day
  starts = Array(ideal_start:Second(cc_step):E)
  ends = starts .+ Second(cc_len) .- Millisecond(convert(Int,1. / fs * 1e3))
  return d2u(starts[findfirst(x -> x >= S, starts)]), d2u(ends[findlast(x -> x <= E,ends)])
end

function slide_ind(startslice::AbstractFloat,endslice::AbstractFloat,fs::AbstractFloat,t::AbstractArray)
  starttime = t[1,2] * 1e-6
  startind = convert(Int,round((startslice - starttime) * fs,digits=4)) + 1
  endind = convert(Int,round((endslice - starttime) * fs,digits=4)) + 1
  return startind,endind
end
