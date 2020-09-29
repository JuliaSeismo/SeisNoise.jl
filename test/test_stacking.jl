# test stacking

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
    Cnew = stack(C)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)

    # test hour stack -> this does nothing
    Cnew = stack(C,interval=Hour(1))
    @test Cnew.t == C.t
    @test Cnew.corr == C.corr

    # test 2hour stack
    Cnew = stack(C,interval=Hour(2))
    @test Cnew.t == collect(0.:7200.:(Nwin-1)*3600.) .+ starttime
    @test size(Cnew.corr) == (N,5)

    # test 2hour stack w/ robuststack
    Cnew = stack(C,interval=Hour(2),stacktype=robuststack)
    @test Cnew.t == collect(0.:7200.:(Nwin-1)*3600.) .+ starttime
    @test size(Cnew.corr) == (N,5)

    # test allstack
    Cnew = stack(C,allstack=true)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)

    # test allstack w/ robuststack
    Cnew = stack(C,allstack=true,stacktype=robuststack)
    @test Cnew.t == [starttime]
    @test size(Cnew.corr) == (N,1)

    # test inplace
    Cnew = stack(C)
    stack!(C)
    @test Cnew.corr == C.corr
    @test Cnew.t == C.t
end

@testset "robust stack" begin
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

@testset "phase weighted stack " begin
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
