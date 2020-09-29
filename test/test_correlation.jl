# test correlation

T = Float32
fs = 100.  # sampling rate
N = 2^14   # number of points
Nwin = 10
cc_len = 40.96 # length of noise window
cc_step = 40.96 # step between windows
freqmin = 10. # minimum frequency
freqmax = 20. # maximum frequency
corners = 4   # number of corners
zerophase = true # do zerophase filtering
maxlag = 20. # seconds of coda to save
starttime = d2u(DateTime(Date(now()))) # starttime for SeisChanel / SeisData
rand_net()    = uppercase(String(rand('A':'Z', 2)))
rand_sta()    = String(rand('A':'Z', rand(3:5)))
rand_loc()    = rand() < 0.5 ? "" : lpad(rand(0:9), 2, "0")
rand_id() = string(rand_net(), ".", rand_sta(), ".", rand_loc(), ".HH", rand('A':'Z'))

@testset "clean up" begin
    A = rand(T,N,Nwin)
    B = clean_up(A,freqmin,freqmax,fs)
    C = clean_up(A,Int(freqmin),Int(freqmax),Int(fs))
    @test size(B) == size(A)
    @test eltype(B) == eltype(A)
    @test eltype(C) == eltype(A)
    @test B != A
    @test C == B

    # test in-place
    clean_up!(A,freqmin,freqmax,fs)
    @test A == B
    @test A == C

    # test RawData
    R = RawData()
    R.x = rand(T,N,Nwin)
    R.fs = fs
    R.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime
    Rclean = clean_up(R,freqmin,freqmax)
    RInt = clean_up(R,Int(freqmin),Int(freqmax))
    @test size(Rclean.x) == size(R.x)
    @test eltype(Rclean.x) == eltype(R.x)
    @test size(RInt.x) == size(R.x)
    @test eltype(RInt.x) == eltype(R.x)
    @test Rclean.x != R.x
    @test Rclean.freqmin == freqmin
    @test Rclean.freqmax == freqmax
    @test RInt.freqmin == freqmin
    @test RInt.freqmax == freqmax

    # test in-place
    clean_up!(R,freqmin,freqmax)
    @test R.x == Rclean.x
    @test R.x == RInt.x
    @test R.freqmin == freqmin
    @test R.freqmax == freqmax

    # test CorrData
    C = CorrData()
    C.corr = rand(T,N,Nwin)
    C.fs = fs
    C.t = collect(0.:3600.:(Nwin-1)*3600.) .+ starttime
    Cclean = clean_up(C,freqmin,freqmax)
    CInt = clean_up(C,Int(freqmin),Int(freqmax))
    @test size(Cclean.corr) == size(C.corr)
    @test eltype(Cclean.corr) == eltype(C.corr)
    @test Cclean.corr != C.corr
    @test Cclean.freqmin == freqmin
    @test Cclean.freqmax == freqmax
    @test CInt.freqmin == freqmin
    @test CInt.freqmax == freqmax 

    # test in-place
    clean_up!(C,freqmin,freqmax)
    @test C.corr == Cclean.corr
    @test C.corr == CInt.corr
    @test C.freqmin == freqmin
    @test C.freqmax == freqmax
end

@testset "correlation" begin
    # test correlate function on arrays
    A = rand(T,N,Nwin) .- T(0.5)
    FFT = rfft(A,1)

    C = correlate(FFT,FFT,N,Int(maxlag * fs))
    @test size(C) == (Int(2 * maxlag * fs) + 1,Nwin)
    @test eltype(C) == eltype(A)

    # test autocorrelation peak occurs at zero lag
    t = -maxlag:1/fs:maxlag
    maxinds = argmax(C,dims=1)
    @test all([t[maxinds[ii][1]] == 0 for ii in length(maxinds)])

    # test phase correlation on array
    A = rand(T,N,Nwin) .- T(0.5)
    FFT = fft(A,1)

    C = phasecorrelate(FFT,FFT,N,Int(maxlag * fs))
    @test size(C) == (Int(2 * maxlag * fs) + 1,Nwin)
    @test eltype(C) == eltype(A)

    # test autocorrelation peak occurs at zero lag
    t = -maxlag:1/fs:maxlag
    maxinds = argmax(C,dims=1)
    @test all([t[maxinds[ii][1]] == 0 for ii in length(maxinds)])

    # test correlation with FFTData
    Ch = SeisIO.RandSeis.randSeisChannel(c=false,s=true)
    ungap!(Ch)
    Ch.x = rand(T,Int(cc_len*Nwin*fs +1)) .- T(0.5)
    Ch.fs = fs
    Ch.t = [1 starttime * 1e6;length(Ch.x) 0]
    Ch.loc = GeoLoc()
    R = RawData(Ch,cc_len,cc_step)
    F = rfft(R)
    C = correlate(F,F,maxlag)
    CInt = correlate(F,F,Int(maxlag))
    @test isa(C,CorrData) # test return type
    @test C.comp == repeat(Ch.id[end],2) # test component name
    @test size(C.corr) == (Int(2 * maxlag * fs) + 1,Nwin) # test size
    @test eltype(C.corr) == eltype(Ch.x) # test type stability
    maxinds = argmax(C.corr,dims=1)
    @test all([t[maxinds[ii][1]] == 0 for ii in length(maxinds)]) # test max args
    @test C.rotated == false # test rotation
    @test C.corr_type == "CC" # test cross-correlation
    @test C.maxlag == CInt.maxlag 
    @test C.corr == CInt.corr

    # test windows that do not overlap
    F1 = deepcopy(F)
    F2 = deepcopy(F)
    F2.t[end] = rand(eltype(F2.t))
    C1 = correlate(F1,F2,maxlag)
    @test size(C1.corr) == (Int(2 * maxlag * fs) + 1,Nwin-1) # test size

    # test with no matching windows
    F2.t = rand(eltype(F2.t),size(F2.fft,2))
    @test_throws ArgumentError correlate(F1,F2,maxlag)

    # test phase correlation
    P = phase(R)
    PCC = correlate(P,P,maxlag,corr_type="PCC")
    @test isa(PCC,CorrData) # test return type
    @test PCC.comp == repeat(Ch.id[end],2) # test component name
    @test size(PCC.corr) == (Int(2 * maxlag * fs) + 1,Nwin) # test size
    @test eltype(PCC.corr) == eltype(Ch.x) # test type stability
    maxinds = argmax(PCC.corr,dims=1)
    @test all([t[maxinds[ii][1]] == 0 for ii in length(maxinds)]) # test max args
    @test PCC.rotated == false # test rotation
    @test PCC.corr_type == "PCC" # test cross-correlation

    # test unrecognized type of correlation
    @test_throws ArgumentError correlate(F,F,maxlag,corr_type="AUTO")

end

@testset "whitening" begin

    pad = 50
    padfreq = pad / fs
    # test whitening array
    A = rand(T,N,Nwin) .- T(0.5)
    F = rfft(A,1)
    Fwhite = whiten(F,freqmin,freqmax,fs,N,pad=pad)

    # test whitening
    freqs = FFTW.rfftfreq(N,fs)
    lowband = findall(0 + eps(T) .< freqs .< freqmin - padfreq)
    passband = findall(freqmin .<= freqs .<= freqmax)
    highband = findall(freqs .> freqmax + padfreq)

    # test filtering
    @test all(Fwhite[lowband,:] .== 0)
    @test all(Fwhite[passband,:] .!= 0)
    @test all(Fwhite[highband,:] .== 0)
    @test size(Fwhite) == size(F)
    @test eltype(Fwhite) == eltype(F)

    # test with low whitening passband
    Flow = whiten(F,0.,freqmax,fs,N,pad=pad)
    @test all(abs.(Flow[1:pad,:]) .< 1.)

    # test with high whitening passband
    Fhigh = whiten(F,freqmin,fs/2,fs,N,pad=pad)
    @test all(abs.(Fhigh[end-pad:end,:]) .< 1.)

    # test with integer arguments 
    FInt = whiten(F,Int(freqmin),Int(freqmax),Int(fs),N,pad=pad)

    # test in-place
    whiten!(F,freqmin,freqmax,fs,N,pad=pad)
    @test F == Fwhite
    @test F == FInt

    # test FFTData
    Ch = SeisIO.RandSeis.randSeisChannel(c=false,s=true)
    Ch.loc = GeoLoc()
    ungap!(Ch)
    Ch.x = rand(T,Int(cc_len*Nwin*fs +1)) .- T(0.5)
    Ch.fs = fs
    Ch.t = [1 starttime * 1e6;length(Ch.x) 0]
    R = RawData(Ch,cc_len,cc_step)
    F = rfft(R)

    Fwhite = whiten(F,freqmin,freqmax)
    @test isa(Fwhite,FFTData)
    @test size(Fwhite.fft) == size(F.fft)
    @test Fwhite.whitened == true
    @test Fwhite.freqmin == freqmin
    @test Fwhite.freqmax == freqmax
    @test Fwhite.fft != F.fft

    # test low/high frequency
    Flowhigh = whiten(F,0.,fs)
    @test Flowhigh.whitened == true
    @test Flowhigh.freqmin == F.freqmin
    @test Flowhigh.freqmax == F.freqmax

    # test low frequency
    Flow = whiten(F,0.,freqmax)
    @test Flow.whitened == true
    @test Flow.freqmin == F.freqmin
    @test Flow.freqmax == freqmax

    # test low frequency
    Fhigh = whiten(F,freqmin,fs)
    @test Fhigh.whitened == true
    @test Fhigh.freqmin == freqmin
    @test Fhigh.freqmax == F.freqmax

    # test with integer input 
    FInt = whiten(F,Int(freqmin),Int(freqmax))
    @test FInt.freqmin == freqmin
    @test FInt.freqmax == freqmax

    # test in-place
    whiten!(F,freqmin,freqmax)
    @test F.fft == Fwhite.fft
    @test F.whitened == true
    @test F.freqmin == freqmin
    @test F.freqmax == freqmax
    @test FInt == F

    # test RawData
    Rwhite = whiten(R,freqmin,freqmax)
    @test isa(Rwhite,RawData)
    @test size(Rwhite.x) == size(R.x)
    @test Rwhite.whitened == true
    @test Rwhite.freqmin == freqmin
    @test Rwhite.freqmax == freqmax
    @test Rwhite.x != R.x

    # test low/high frequency
    Rlowhigh = whiten(R,0.,fs)
    @test Rlowhigh.whitened == true
    @test Rlowhigh.freqmin == R.freqmin
    @test Rlowhigh.freqmax == R.freqmax

    # test low frequency
    Rlow = whiten(R,0.,freqmax)
    @test Rlow.whitened == true
    @test Rlow.freqmin == R.freqmin
    @test Rlow.freqmax == freqmax

    # test low frequency
    Rhigh = whiten(R,freqmin,fs)
    @test Rhigh.whitened == true
    @test Rhigh.freqmin == freqmin
    @test Rhigh.freqmax == R.freqmax

    # test with integer input 
    RInt = whiten(R,Int(freqmin),Int(freqmax))
    @test RInt.freqmin == freqmin 
    @test RInt.freqmax == freqmax

    # test in-place
    whiten!(R,freqmin,freqmax)
    @test R.x == Rwhite.x
    @test R.whitened == true
    @test R.freqmin == freqmin
    @test R.freqmax == freqmax
    @test RInt == R

    # test coherence with smoothing
    R = RawData(Ch,cc_len,cc_step)
    F = rfft(R)
    Fcoh = coherence(F,20,0.01)
    @test isa(Fcoh,FFTData)
    @test size(Fcoh.fft) == size(F.fft)
    @test eltype(Fcoh.fft) == eltype(F.fft)
    @test Fcoh.fft != F.fft

    # test coherence no smoothing
    Fcoh = coherence(F,20)
    @test isa(Fcoh,FFTData)
    @test size(Fcoh.fft) == size(F.fft)
    @test eltype(Fcoh.fft) == eltype(F.fft)
    @test Fcoh.fft != F.fft

    # test inplace
    coherence!(F,20)
    @test F.fft == Fcoh.fft

    # test deconvolution with smoothing
    R = RawData(Ch,cc_len,cc_step)
    F = rfft(R)
    Fdec = deconvolution(F,20,0.01)
    @test isa(Fdec,FFTData)
    @test size(Fdec.fft) == size(F.fft)
    @test eltype(Fdec.fft) == eltype(F.fft)
    @test Fdec.fft != F.fft

    # test deconvolution no smoothing
    Fdec = deconvolution(F,20)
    @test isa(Fdec,FFTData)
    @test size(Fdec.fft) == size(F.fft)
    @test eltype(Fdec.fft) == eltype(F.fft)
    @test Fdec.fft != F.fft

    # test inplace
    deconvolution!(F,20)
    @test F.fft == Fdec.fft
end

@testset "corr mapping" begin
    Ch = SeisIO.RandSeis.randSeisChannel(c=false,s=true)
    ungap!(Ch)
    Ch.x = rand(T,Int(cc_len*Nwin*fs +1)) .- T(0.5)
    Ch.fs = fs
    Ch.t = [1 starttime * 1e6;length(Ch.x) 0]
    Ch.loc = GeoLoc()
    R = RawData(Ch,cc_len,cc_step)
    Fs = [rfft(R) for ii in 1:5]

    # change FFTData names
    for ii = 1:5
        Fs[ii].name = rand_id()
    end

    # test corrmap with defaults
    tmpdir = tempdir()
    corrmap(deepcopy(Fs),maxlag,tmpdir)
    files = glob("*.jld2",tmpdir)
    @test length(files) == 10
    rm.(files)

    # test with coherence
    corrmap(deepcopy(Fs),maxlag,tmpdir,smooth_type="coherence")
    files = glob("*.jld2",tmpdir)
    @test length(files) == 10
    rm.(files)

    # test with deconvolution
    corrmap(deepcopy(Fs),maxlag,tmpdir,smooth_type="deconvolution")
    files = glob("*.jld2",tmpdir)
    @test length(files) == 10
    rm.(files)

    # test with stacking
    corrmap(deepcopy(Fs),maxlag,tmpdir,interval=Day(1))
    files = glob("*.jld2",tmpdir)
    @test length(files) == 10
    rm.(files)
end
