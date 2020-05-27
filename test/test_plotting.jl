# test plotting

@testset "corr plot" begin
    T = Float32
    fs = 100.  # sampling rate
    Nwin = 10
    maxlag = 20. # seconds of coda to save
    starttime = d2u(DateTime(Date(now()))) # starttime for SeisChanel / SeisData
    C = CorrData()
    t = T.(-maxlag:1/fs:maxlag)
    C.corr = rand(T,length(t),Nwin) .+ T(5) .* repeat(sin.(t),inner=(1,Nwin))
    C.fs = fs
    C.maxlag = maxlag
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime
    @test !isa(try corrplot(C) catch ex ex end, Exception)    # Passes
end
