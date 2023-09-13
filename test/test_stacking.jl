# testSeisNoise.stacking

T = Float32
fs = 100.  # sampling rate
N = 2^14   # number of points
Nwin = 10
freqmin = 10. # minimum frequency
freqmax = 20. # maximum frequency
corners = 4   # number of corners
zerophase = true # do zerophase filtering
starttime = d2u(DateTime(Date(now()))) # starttime for SeisChanel / SeisData

@testset "stack" begin
    A = rand(T,N,Nwin)
    C = CorrData()
    C.corr = A
    C.fs = fs
    # corrleation every hour
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime

    # test daystack
    Cnew = SeisNoise.stack(C)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)

    # test hourSeisNoise.stack -> this does nothing
    Cnew = SeisNoise.stack(C,interval=Hour(1))
    @test Cnew.t == C.t
    @test Cnew.corr == C.corr

    # test 2hourSeisNoise.stack
    Cnew = SeisNoise.stack(C,interval=Hour(2))
    @test Cnew.t == collect(0.:7200.:(Nwin-1)*3600.) .+ starttime
    @test size(Cnew.corr) == (N,5)

    # test 2hourSeisNoise.stack w/ robuststack
    Cnew = SeisNoise.stack(C,interval=Hour(2),stacktype=robuststack)
    @test Cnew.t == collect(0.:7200.:(Nwin-1)*3600.) .+ starttime
    @test size(Cnew.corr) == (N,5)

    # test allstack
    Cnew = SeisNoise.stack(C,allstack=true)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)

    # test allstack w/ robuststack
    Cnew = SeisNoise.stack(C,allstack=true,stacktype=robuststack)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)

    # test inplace
    Cnew = SeisNoise.stack(C)
   SeisNoise.stack!(C)
    @test Cnew.corr == C.corr
    @test Cnew.t == C.t
end

@testset "robustSeisNoise.stack" begin
    A = rand(T,N,Nwin) .+ T.(sin.(range(0.,10π,length=N)))
    B = robuststack(A)
    @test size(B) == (N,1)
    @test eltype(B) == eltype(A)

    # test convergence
    B = robuststack(A,ϵ=Float32(1e-3))
    @test size(B) == (N,1)
    @test eltype(B) == eltype(A)

    # test maxiter
    B = robuststack(A,maxiter=15)
    @test size(B) == (N,1)
    @test eltype(B) == eltype(A)

    # test CorrData
    C = CorrData()
    C.corr = A
    C.fs = fs
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime
    Cnew = robuststack(C)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)
    @test eltype(Cnew.corr) == eltype(A)

    # test in-place
    robuststack!(C)
    @test isapprox(Cnew.corr,C.corr)
    @test Cnew.t == C.t
end

@testset "median mute" begin
    A = rand(T,N,Nwin)
    high = 2.
    low = 0.5
    # set two out of 10 to have large amplute
    A[:,1] .*= 1000
    # set two out of 10 to have small amplitude
    A[:,2] ./= 1000

    # test indexes for 0 < medianamp < high
    ind = SeisNoise.medianmuteind(A, high)
    @test length(ind) == Nwin - 1

    # test indexes for low < medianamp < high
    ind = SeisNoise.medianmuteind(A, high, low)
    @test length(ind) == Nwin - 2

    # test when low > high
    @test_throws AssertionError SeisNoise.medianmuteind(A, low, high)

    # test with Nans
    # set two columns to have NaNs
    A[rand(1:N,3),3] .= NaN
    ind = SeisNoise.medianmuteind(A, high, low)
    @test length(ind) == Nwin - 3

    # test with Infs
    # set two columns to have Infs
    A[rand(1:N,3),4] .= Inf
    ind = SeisNoise.medianmuteind(A, high, low)
    @test length(ind) == Nwin - 4

    # test with CorrData
    C = CorrData()
    C.corr = A
    C.fs = fs
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime

    # test allocating version
    Cnew = medianmute(deepcopy(C), high)
    @test Cnew != C
    @test size(Cnew.corr) == (N,7)
    @test length(Cnew.t) == 7

    # test with high and low thresholds
    Cnew = medianmute(deepcopy(C), high, low)
    @test Cnew != C
    @test size(Cnew.corr) == (N,6)
    @test length(Cnew.t) == 6

    # test in-place version
    medianmute!(C, high, low)
    @test C == Cnew
end

@testset "phase weightedSeisNoise.stack " begin
    A = rand(T,N,Nwin) .+ T.(sin.(range(0.,10π,length=N)))
    B = pws(A)
    @test size(B) == (N,1)
    @test eltype(B) == eltype(A)

    # test power
    B = pws(A,pow=1)
    @test size(B) == (N,1)
    @test eltype(B) == eltype(A)

    # test CorrData
    C = CorrData()
    C.corr = A
    C.fs = fs
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime
    Cnew = pws(C)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)
    @test eltype(Cnew.corr) == eltype(A)

    # test in-place
    pws!(C,pow=2.)
    @test isapprox(Cnew.corr,C.corr)
    @test Cnew.t == C.t
end

@testset "remove nan" begin
    C = CorrData()
    A = rand(T,N,Nwin)
    C.corr = A
    C.fs = fs
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime

    # add nan to 4 columns
    C.corr[unique(rand(1:N,100)),1:3:end] .= NaN

    # remove nan
    Cnew = remove_nan(C)
    @test size(Cnew.corr) == (N,6)
    @test eltype(Cnew.corr) == eltype(A)
    @test length(Cnew.t) == 6

    # test inplace
    remove_nan!(C)
    @test Cnew.corr == C.corr
    @test Cnew.t == C.t

    # test no-nan
    C = CorrData()
    A = rand(T,N,Nwin)
    C.corr = A
    C.fs = fs
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime
    Cnew = remove_nan(C)
    @test Cnew.corr == C.corr
    @test Cnew.t == C.t

    # test all NaN
    C.corr[:,:] .= NaN
    @test_throws ArgumentError remove_nan(C)
end

@testset "smooth CorrData" begin
    C = CorrData()
    A = rand(T,10001,Nwin)
    C.corr = A
    # time vector is  [T00:00:00, T01:00:00, T02:00:00, ...]
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime

    # smooth every 2 hours
    Cnew = smooth(C,Hour(2))
    # test that 2nd correlation is mean of first two hours
    @test Cnew.corr[:,2] == vec(mean(C.corr[:,1:2],dims=2))

    # add two hours to last 5 correlations
    C.t[3:10] .+= 3600 * 3
    Cnew = smooth(C,Hour(3))
    # test that 2nd correlation is mean of first two hours
    @test Cnew.corr[:,2] == vec(mean(C.corr[:,1:2],dims=2))
    # test that 3rd correlation remains unchanged
    @test Cnew.corr[:,3] == C.corr[:,3]
    # test that 5th correlation is mean of 3:5
    @test Cnew.corr[:,5] == vec(mean(C.corr[:,3:5],dims=2))

    # randomize starttimes
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime
    ind = SeisNoise.sample(1:Nwin,Nwin,replace=false)
    Cnew = smooth(C[ind],Hour(2))

    # test correlations are ordered
    @test all(diff(Cnew.t) .== 3600)
    @test Cnew.corr[:,2] == vec(mean(C.corr[:,1:2],dims=2))

    # change time vector to Days
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .* 24 .+ starttime
    Cnew = smooth(C,Day(3))
    @test Cnew.corr[:,3] == vec(mean(C.corr[:,1:3],dims=2))

    # test inplace version
    smooth!(C,Day(3))
    @test Cnew == C
end

@testset "shorten CorrData" begin
    C = CorrData()
    A = rand(T,10001,Nwin)
    C.corr = A
    C.fs = fs
    C.maxlag = 50.
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime

    # test shorten to 40s
    newlag = 40.
    Cnew = shorten(C,newlag)
    @test size(Cnew.corr) == (Int(newlag * fs * 2 + 1),Nwin)
    @test eltype(Cnew.corr) == eltype(C.corr)
    @test Cnew.maxlag == newlag

    # test with Int
    CInt = shorten(C,Int(newlag))
    @test size(CInt.corr) == (Int(newlag * fs * 2 + 1),Nwin)
    @test eltype(CInt.corr) == eltype(C.corr)
    @test CInt.maxlag == newlag

    # test newlag > maxlag
    newlag = 100.
    @test_throws AssertionError shorten(C,newlag)

    # test newlag =  0
    newlag = 0.
    @test_throws AssertionError shorten(C,newlag)

    # test newlag <=  0
    newlag = -5.
    @test_throws AssertionError shorten(C,newlag)

    # test where newlag % (1/C.fs) != 0
    newlag = π
    Cpi = shorten(C,π)
    @test Cpi.maxlag == 3.14

    # test in-place
    newlag = 40.
    Cnew = shorten(C,newlag)
    shorten!(C,newlag)
    @test Cnew.corr == C.corr
    @test Cnew.maxlag == C.maxlag
end

@testset "Adaptive Filter" begin
    # test adaptive filter
    A = rand(T,N,Nwin)
    B = deepcopy(A)
    window=10.
    adaptive_filter!(B,window,fs)
    @test size(B) == size(A)
    @test eltype(B) == eltype(B)

    # test overlap
    B = deepcopy(A)
    @test_throws ArgumentError adaptive_filter!(B,window,fs,overlap=0.)

    @test_throws ArgumentError adaptive_filter!(B,window,fs,overlap=1.)

    # test 1D vector
    v = rand(T,N)
    vcopy = deepcopy(v)
    adaptive_filter!(vcopy,window,fs)
    @test v == vcopy

    # test CorrData
    C = CorrData()
    A = rand(T,N,Nwin)
    C.corr = A
    C.fs = fs
    Cnew = adaptive_filter(C,window)
    @test size(Cnew.corr) == size(C.corr)
    @test eltype(Cnew.corr) == eltype(C.corr)

    # test CorrData in-place
    adaptive_filter!(C,window)
    @test Cnew.corr == C.corr

    # test ACF_kernel
    A = rand(T,N,Nwin)
    B = SeisNoise.ACF_kernel(A)
    @test size(B) == size(A)
    @test eltype(B) == eltype(A)

end

@testset "robustpws" begin
    A = rand(T,N,Nwin) .+ T.(sin.(range(0.,10π,length=N)))
    B = robustpws(A)
    @test size(B) == (N,1)
    @test eltype(B) == eltype(A)

    # test CorrData
    C = CorrData()
    C.corr = A
    C.fs = fs
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime
    Cnew = robustpws(C)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)
    @test eltype(Cnew.corr) == eltype(A)

    # test in-place
    robustpws!(C)
    @test Cnew.corr == C.corr
    @test Cnew.t == C.t
end
