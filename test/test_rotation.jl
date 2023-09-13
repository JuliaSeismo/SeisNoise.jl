# test rotatation tools

T = Float32
fs = 100.  # sampling rate
N = 2^14   # number of points
Nwin = 10
cc_len = 40.96 # length of noise window
cc_step = 40.96 # step between windows
azi = 320 / 180 * π
baz = 140 / 180 * π
maxlag = 20. # seconds of coda to save
starttime = d2u(DateTime(Date(now()))) # starttime for SeisChanel / SeisData

@testset "rotation" begin

    Ch = SeisBase.RandSeis.randSeisChannel(c=false,s=true)
    ungap!(Ch)
    Ch.x = rand(T,Int(cc_len*Nwin*fs +1)) .- T(0.5)
    Ch.fs = fs
    Ch.t = [1 starttime * 1e6;length(Ch.x) 0]
    Ch.loc = GeoLoc()
    R = RawData(Ch,cc_len,cc_step)
    F1s = [rfft(R) for ii in 1:3]
    F2s = [rfft(R) for ii in 1:3]

    # change FFTData names to ENZ naming convention
    F1s[1].name = "NN.STA1.00.HHE"
    F1s[2].name = "NN.STA1.00.HHN"
    F1s[3].name = "NN.STA1.00.HHZ"
    F2s[1].name = "NN.STA2.00.HHE"
    F2s[2].name = "NN.STA2.00.HHN"
    F2s[3].name = "NN.STA2.00.HHZ"

    # cross-correlate
    Cs = Array{CorrData}(undef,9)
    let  kk = 1
        for ii = 1:3
            for jj = 1:3
                Cs[kk] = correlate(F1s[ii],F2s[jj],maxlag)
                kk += 1
            end
        end
    end

    # rotate
    Crot = rotate(Cs,azi,baz)
    comps = [Crot[ii].comp for ii = 1:9]
    @test Set(comps) == Set(["RR","RT","RZ","TR","TT","TZ","ZR","ZT","ZZ"])
    @test all([Crot[ii].rotated == true for ii = 1:9])
    @test all([eltype(Crot[ii].corr) == T for ii = 1:9])

    # NEEDED: test with real data

    # test in-place
    rotate!(Cs,azi,baz)
    comps = [Cs[ii].comp for ii = 1:9]
    @test Set(comps) == Set(["RR","RT","RZ","TR","TT","TZ","ZR","ZT","ZZ"])
    @test all([Cs[ii].rotated == true for ii = 1:9])
    @test all([Cs[ii].corr == Crot[ii].corr for ii = 1:9])
    @test all([eltype(Cs[ii].corr) == T for ii = 1:9])
end
