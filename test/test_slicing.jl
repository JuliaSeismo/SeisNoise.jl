# test slicing of Data

C = SeisIO.RandSeis.randSeisChannel(s=true,c=false)
S = SeisIO.RandSeis.randSeisData(s=1.,c=0.)[1:1]
ungap!(C)
ungap!(S)

@testset "start_end" begin
@test eltype(start_end(C)) == DateTime
@test isa(start_end(C),Array{DateTime})
@test size(start_end(C),1) + 1 == size(C.t,1)
@test eltype(start_end(S)) == Array{DateTime,1}
# start_end(S) might need a revamp
end

@testset "slide" begin
# 1000 second long signal
A = rand(10000)
fs = 10.
cc_step = cc_len = 100.
starttime = d2u(DateTime(Date(now())))
endtime = starttime + length(A) / fs - 1 / fs

# test even splitting
@test slide(A,cc_len,cc_step,fs,starttime,endtime)[1] == reshape(A,Int(fs*cc_len),Int(length(A) / cc_len / fs) )
@test slide(A,cc_len,cc_step,fs,starttime,endtime)[2] == Array(range(starttime,stop=endtime,step=cc_step))

# test uneven splitting
cc_len = cc_step = 33.
num_windows = Int(length(A) ÷ cc_len ÷ fs)
pts = Int(cc_len * fs * num_windows )
@test slide(A,cc_len,cc_step,fs,starttime,endtime)[1] == reshape(A[1:pts],Int(cc_len * fs),num_windows)
@test slide(A,cc_len,cc_step,fs,starttime,endtime)[2] == Array(range(starttime,stop=endtime,step=cc_step))[1:num_windows]

# test overlapping windows
cc_len = 100.
cc_step = 50.
out, starts = slide(A,cc_len,cc_step,fs,starttime,endtime)
# test first column
@test out[:,1] == A[1:Int(cc_len*fs)]
# test overlap
@test out[Int(cc_step*fs)+1:end,1] == out[1:Int(cc_step*fs),2]
@test length(starts) == size(out,2)
@test eltype(out) == eltype(A)

# test SeisChannel slide
C = SeisChannel()
C.x = A
C.fs = fs
C.t = [1 starttime * 1e6;length(C.x) 0]
@test slide(C,cc_len,cc_step)[1] == out
end

# test nearest_start_end function
@testset "nearest_start_end" begin
A = rand(10000)
fs = 10.
# start channel at 00:01:12
td = DateTime(Date(now()))
starttime = td + Minute(1) +  Second(12)
endtime = starttime + Second(length(A) / fs)
C = SeisChannel()
C.x = A
C.fs = fs
C.t = [1 d2u(starttime) * 1e6;length(C.x) 0]
cc_len = cc_step = 100.
# with 100s step, nearest start is 00:01:40
# nearest_end is 1000s from 00:00:00 = 00:16:40 - 1 / fs = 00:16:39
nextstart = d2u(td + Minute(1) +  Second(40))
lastend = d2u(td + Minute(16) + Second(39) + Millisecond(900))
@test nearest_start_end(C,cc_len,cc_step) == (nextstart,lastend)
@test eltype(nearest_start_end(C,cc_len,cc_step)) == Float64

# test with DateTime creation
@test nearest_start_end(starttime,endtime,fs,cc_len,cc_step) == (nextstart,lastend)
@test nearest_start_end(starttime,endtime,fs,cc_len,cc_step) == nearest_start_end(C,cc_len,cc_step)
end

@testset "slide_ind" begin
A = rand(10000)
fs = 10.
# start at midnight
td = DateTime(Date(now()))
C = SeisChannel()
C.x = A
C.fs = fs
C.t = [1 d2u(td) * 1e6;length(C.x) 0]
cc_len = cc_step = 100.
startslice, endslice = nearest_start_end(C,cc_len,cc_step)
@test slide_ind(startslice,endslice,fs,C.t) == (1,10000)

# start channel at 00:01:12
starttime = td + Minute(1) +  Second(12)
endtime = starttime + Second(length(A) / fs) - Millisecond(1000 * 1 / fs)
C.t = [1 d2u(starttime) * 1e6;length(C.x) 0]
startslice, endslice = nearest_start_end(C,cc_len,cc_step)
# start at 00:01:12
# want to start at 00:01:40
# or 28 seconds + 1 sample == 281 samples
# endtime is 00:16:39.9 or 15:27.9 or 9280 samples
@test slide_ind(startslice,endslice,fs,C.t) == (281,9280)
end
