# test deprecated functions

T = Float32
fs = 100.  # sampling rate
N = 2^14   # number of points
Nwin = 10
freqmin = 10. # minimum frequency
freqmax = 20. # maximum frequency
corners = 4   # number of corners
zerophase = true # do zerophase filtering
maxlag = 20.     # maximum lag time in coda
starttime = d2u(DateTime(Date(now()))) # starttime for SeisChanel / SeisData


# test FFTData creation with compute_fft
@testset "compute fft" begin

    cc_len = 10.24
    cc_step = 5.12
    C = SeisChannel()
    C.x = rand(T,N)
    C.fs = fs
    C.id = "MZ.KTFD.08.ENZ"
    C.t = [1 (starttime + 0.004) * 1e6;N 0]
    R = RawData(C,cc_len,cc_step)
    R.gain = 1.234234e6
    F = compute_fft(R)
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

# test correlation with FFTData
@testset "compute cc" begin
    cc_len = 40.96 # length of noise window
    cc_step = 40.96 # step between windows
    Ch = SeisIO.RandSeis.randSeisChannel(c=false,s=true)
    ungap!(Ch)
    Ch.x = rand(T,Int(cc_len*Nwin*fs +1)) .- T(0.5)
    Ch.fs = fs
    Ch.t = [1 starttime * 1e6;length(Ch.x) 0]
    Ch.loc = GeoLoc()
    R = RawData(Ch,cc_len,cc_step)
    F = rfft(R)
    C = compute_cc(F,F,maxlag)
    @test isa(C,CorrData) # test return type
    @test C.comp == repeat(Ch.id[end],2) # test component name
    @test size(C.corr) == (Int(2 * maxlag * fs) + 1,Nwin) # test size
    @test eltype(C.corr) == eltype(Ch.x) # test type stability
    t = -maxlag:1/fs:maxlag
    maxinds = argmax(C.corr,dims=1)
    @test all([t[maxinds[ii][1]] == 0 for ii in length(maxinds)]) # test max args
    @test C.rotated == false # test rotation
    @test C.corr_type == "CC" # test cross-correlation

    # test with no intersecting windows
    F1 = deepcopy(F)
    F2 = deepcopy(F)
    F2.t = rand(eltype(F2.t),size(F2.fft,2))
    @test_throws ArgumentError compute_cc(F1,F2,maxlag)

end
