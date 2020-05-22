# test compute fft

T = Float32
fs = 100.  # sampling rate
N = 2^14   # number of points
Nwin = 10
freqmin = 10. # minimum frequency
freqmax = 20. # maximum frequency
corners = 4   # number of corners
zerophase = true # do zerophase filtering
starttime = d2u(DateTime(Date(now()))) # starttime for SeisChanel / SeisData

# test process raw
@testset "process raw" begin

    S = SeisData(2)
    S.x[1] = rand(T,N)
    S.id[1] = "EO.SYT.02.LLE"
    S.fs[1] = fs
    S.t[1] = [1 starttime * 1e6;N 0]
    S.x[2] = rand(T,N)
    S.id[2] = "MZ.KTFD.08.ENZ"
    S.fs[2] = fs
    S.t[2] = [1 (starttime + 0.004) * 1e6;N 0]

    # test without downsampling
    Snew = process_raw(S,fs)
    @test size(Snew.x[1]) == size(S.x[1])
    @test eltype(Snew.x[1]) == eltype(S.x[1])
    @test size(Snew.x[2]) == size(S.x[2])
    @test eltype(Snew.x[2]) == eltype(S.x[2])
    @test Snew.fs == S.fs
    @test Snew.t[1] == S.t[1]
    @test Snew.t[2][1,2] == S.t[2][1,2] - 0.004 * 1e6 # test phase shift

    # test with downsampling
    Snew = process_raw(S,20.)
    @test size(Snew.x[1]) != size(S.x[1])
    @test eltype(Snew.x[1]) == eltype(S.x[1])
    @test size(Snew.x[2]) != size(S.x[2])
    @test eltype(Snew.x[2]) == eltype(S.x[2])
    @test Snew.fs == [20.,20.]
    @test Snew.t[1][2,1] == length(Snew.x[1])
    @test Snew.t[2][2,1] == length(Snew.x[2])
    @test Snew.t[1] != S.t[1]
    @test Snew.t[2][1,2] == S.t[2][1,2] - 0.004 * 1e6 # test phase shift

    # test without ϕshift
    Sphi = process_raw(S,20.,ϕshift=false)
    @test Sphi.x[1] == Snew.x[1]
    @test Sphi.x[2] != Snew.x[2]
    @test Sphi.t[2] == Snew.t[2]

    # test in-place
    process_raw!(S,20.)
    @test Snew == S

    # test with SeisChannel
    C = SeisChannel()
    C.x = rand(T,N)
    C.fs = fs
    C.t = [1 (starttime + 0.004) * 1e6;N 0]

    # test no downsample
    Cnew = process_raw(C,fs)
    @test size(Cnew.x) == size(C.x)
    @test eltype(Cnew.x) == eltype(C.x)
    @test Cnew.fs == C.fs
    @test Cnew.t[1,2] == C.t[1,2] - 0.004 * 1e6 # test phase shift

    # test downsample
    Cnew = process_raw(C,20.)
    @test size(Cnew.x) != size(C.x)
    @test eltype(Cnew.x) == eltype(C.x)
    @test Cnew.fs == 20.
    @test Cnew.t[1,2] == C.t[1,2] - 0.004 * 1e6 # test phase shift

    # test in-place
    process_raw!(C,20.)
    @test Cnew == C
end

# test FFTData creation with rfft
@testset "rfft" begin

    cc_len = 10.24
    cc_step = 5.12
    C = SeisChannel()
    C.x = rand(T,N)
    C.fs = fs
    C.id = "MZ.KTFD.08.ENZ"
    C.t = [1 (starttime + 0.004) * 1e6;N 0]
    R = RawData(C,cc_len,cc_step)
    R.gain = 1.234234e6
    F = rfft(R)
    @test F.name == C.id == R.name
    @test F.id == R.id
    @test F.loc == R.loc
    @test F.fs == R.fs
    @test F.gain == R.gain
    @test F.freqmin == R.freqmin
    @test F.freqmax == F.freqmax
    @test F.cc_len == R.cc_len == cc_len
    @test F.cc_step == R.cc_step == cc_step
    @test F.whitened == R.whitened
    @test F.time_norm == R.time_norm
    @test F.resp == R.resp
    @test F.t == R.t
    @test length(F.t) == size(F.fft,2)
end

@testset "phase" begin

    # test with all zero row
    A = rand(T,N,Nwin)
    A[:,1] .= T(0)
    p = phase(A)
    @test all(.!SeisNoise.isweird.(p)) # test no NaNs

    # test phase
    A = rand(T,N,Nwin)
    p = phase(A)
    @test size(A) == size(p)
    @test all(p[N ÷ 2 + 2 + isodd(N),1] .== T(0))

    # test phase RawData
    R = RawData()
    R.x = A
    R.fs = fs
    F = phase(R)
    @test isa(F,FFTData)
    @test all(.!SeisNoise.isweird.(F.fft))
    @test F.fft == p
end

@testset "mute" begin

    # generate random signal with zero mean
    A = rand(T,N,Nwin)
    A .-= mean(A,dims=1)
    A[1000:1200,:] .*= 1000
    M = mute(A)
    @test size(M) == size(A)
    @test eltype(M) == eltype(A)

    # test begin , middle, end changes
    @test all(M[1:500,:] .== A[1:500,:])
    @test median(M[1000:1200,:]) == T(0)
    @test all(M[end-500:end,:] .== A[end-500:end,:])

    # test in-place
    mute!(A)
    @test A == M

    # test RawData
    R = RawData()
    A = rand(T,N,Nwin)
    A .-= mean(A,dims=1)
    A[1000:1200,:] .*= 1000
    R.x = A
    M = mute(R)
    @test isa(M,RawData)
    @test size(M.x) == size(A)
    @test eltype(M.x) == eltype(A)

    # test in-place
    mute!(R)
    @test R.x == M.x
end

@testset "clip/clamp" begin
    Nspikes = 100
    A = rand(T,N,Nwin) .- 0.5
    ind = rand(1:N,Nspikes)
    A[ind,:] .= randn(Nspikes,Nwin) .* 100000
    A[:,1] .= T(0)

    # test clipping
    B = clip(A,3)
    @test size(B) == size(A)
    @test eltype(B) == eltype(A)
    @test B != A
    @test all(maximum(abs.(B),dims=1) .<= maximum(abs.(A),dims=1))

    # test in-place
    clip!(A,3)
    @test A == B

    # test 1D array
    A = rand(T,N) .- 0.5
    ind = rand(1:N,Nspikes)
    A[ind] .= randn(Nspikes) .* 100000
    B = clip(A,3)
    @test size(B) == size(A)
    @test eltype(B) == eltype(A)
    @test all(maximum(abs.(B),dims=1) .<= maximum(abs.(A),dims=1))

    # test RawData
    Nspikes = 100
    A = rand(T,N,Nwin) .- 0.5
    ind = rand(1:N,Nspikes)
    A[ind,:] .= randn(Nspikes,Nwin) .* 100000
    A[:,1] .= T(0)
    R = RawData()
    R.x = A
    Rclip = clip(R,3)
    @test isa(Rclip,RawData)
    @test size(Rclip.x) == size(R.x)
    @test eltype(Rclip.x) == eltype(R.x)
    @test Rclip.x != R.x
    @test all(maximum(abs.(Rclip.x),dims=1) .<= maximum(abs.(R.x),dims=1))

    # test in-place
    clip!(R,3)
    @test R.x == Rclip.x

    # test clamping w RawData
    A = rand(T,N,Nwin) .- 0.5
    ind = rand(1:N,Nspikes)
    A[ind,:] .= randn(Nspikes,Nwin) .* 100000
    A[:,1] .= T(0)
    R = RawData()
    R.x = deepcopy(A)
    Rclamp = clamp(R,π)
    @test isa(Rclamp,RawData)
    @test size(Rclamp.x) == size(R.x)
    @test eltype(Rclamp.x) == eltype(R.x)
    @test all(abs.(Rclamp.x) .<= π)
    @test Rclamp.x != R.x

    # test in place
    clamp!(R,π)
    @test R.x == Rclamp.x

    # test clamping with two values
    R.x = deepcopy(A)
    Rclamp = clamp(R,-0.75,0.5)
    @test isa(Rclamp,RawData)
    @test size(Rclamp.x) == size(R.x)
    @test eltype(Rclamp.x) == eltype(R.x)
    @test all(Rclamp.x .>= -0.75)
    @test all(Rclamp.x .<= 0.5)

    # test switching hi and lo
    @test_throws ArgumentError clamp(R,0.5,-0.75)

    # test in-place
    clamp!(R,-0.75,0.5)
    @test R.x == Rclamp.x
end

@testset "remove_amp" begin
    # create array with std == 1
    A = rand(T,N,Nwin) .- 0.5
    A ./= std(A,dims=1)

    # add amplitude spikes
    Nspikes = 100
    cols = unique(rand(1:Nwin,Nwin÷2))
    ind = rand(1:N,Nspikes)
    R = RawData()
    R.x = deepcopy(A)
    R.x[ind,cols] .= randn(Nspikes,length(cols)) .* 100000
    R.t = collect(1.:10.)
    Ramp = remove_amp(R)
    @test isa(Ramp,RawData)
    @test size(Ramp.x,2) == Nwin - length(cols)
    @test length(Ramp.t) == Nwin - length(cols)

    # test in-place
    remove_amp!(R)
    @test R.x == Ramp.x

    # add amplitude spikes to all columns
    R.x = deepcopy(A)
    R.x[ind,:] .= randn(Nspikes) .* 100000
    @test_throws ErrorException remove_amp(R)

    # test with zeros
    cols = unique(rand(1:Nwin,Nwin÷2))
    R.x = deepcopy(A)
    R.x[:,cols] .= T(0)
    R.t = collect(1.:10.)
    Ramp = remove_amp(R)
    @test isa(Ramp,RawData)
    @test size(Ramp.x,2) == Nwin - length(cols)
    @test length(Ramp.t) == Nwin - length(cols)

    # test with all zeros
    R.x .= zero(T)
    @test_throws ErrorException remove_amp(R)
end

@testset "one-bit" begin
    A = rand(T,N,Nwin) .- T(0.5)
    R = RawData()
    R.x = A

    # one-bit RawData
    R1 = onebit(R)
    @test isa(R1,RawData)
    @test maximum(R1.x) == T(1)
    @test size(R1.x) == size(A)
    @test eltype(R1.x) == eltype(A)
    @test minimum(R1.x) == T(-1)
    @test maximum(R1.x) == T(1)

    # test in-place
    onebit!(R)
    @test R.x == R1.x

end

@testset "nonzero" begin
    A = rand(T,N,Nwin)
    cols = unique(rand(1:Nwin,Nwin÷2))
    A[:,cols] .= T(0)
    ind = SeisNoise.nonzero(A)
    # test no intersect between cols and ind
    @test setdiff(cols,ind) == cols
    # test ind are Int
    @test all(isa.(ind,Int))
end
