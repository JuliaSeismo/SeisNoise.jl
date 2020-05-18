# test filtering

T = Float32
fs = 100.  # sampling rate
N = 2^14   # number of points
Nwin = 10
freqmin = 10. # minimum frequency
freqmax = 20. # maximum frequency
corners = 4   # number of corners
zerophase = true # do zerophase filtering
starttime = d2u(DateTime(Date(now()))) # starttime for SeisChanel / SeisData

# test bandpass
@testset "bandpass" begin

# create white noise
A = rand(T,N,Nwin)
B = bandpass(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
# check size and type stability
@test size(B) == size(A)
@test eltype(B) == eltype(A)
bandpass!(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
@test A == B

# test power in pass bands
freqs = FFTW.rfftfreq(N,fs)
lowband = findall(0 + eps(T) .< freqs .< freqmin)
passband = findall(freqmin .<= freqs .<= freqmax)
highband = findall(freqs .> freqmax)
F = rfft(A,1)
Fpow = abs.(F)
# test lowband
@test all(median(Fpow[lowband,:],dims=1) .< median(Fpow[passband,:],dims=1))

# test highband
@test all(median(Fpow[highband,:],dims=1) .< median(Fpow[passband,:],dims=1))

# test zerophase
A = rand(T,N,Nwin)
@test bandpass(A,freqmin,freqmax,fs,corners=corners,zerophase=true) != bandpass(A,freqmin,
               freqmax,fs,corners=corners,zerophase=false)

# test corners
@test bandpass(A,freqmin,freqmax,fs,corners=4,zerophase=true) != bandpass(A,freqmin,
               freqmax,fs,corners=3,zerophase=true)

# test RawData filtering
R = RawData()
R.x = A
R.fs = fs
Rnew = bandpass(R,freqmin,freqmax)
@test Rnew.fs == R.fs
@test Rnew.freqmin == freqmin
@test Rnew.freqmax == freqmax
bandpass!(R,freqmin,freqmax)
@test R.fs == R.fs
@test R.freqmin == freqmin
@test R.freqmax == freqmax
@test R.x == Rnew.x

# test zerophase
A = rand(T,N,Nwin)
R.x = A
@test bandpass(R,freqmin,freqmax,corners=corners,zerophase=true).x != bandpass(R,freqmin,
               freqmax,corners=corners,zerophase=false).x

# test corners
@test bandpass(R,freqmin,freqmax,corners=4,zerophase=true).x != bandpass(R,freqmin,
              freqmax,corners=3,zerophase=true).x

# test CorrData filtering
C = CorrData()
C.corr = A
C.fs = fs
Cnew = bandpass(C,freqmin,freqmax)
@test Cnew.fs == C.fs
@test Cnew.freqmin == freqmin
@test Cnew.freqmax == freqmax
bandpass!(C,freqmin,freqmax)
@test C.fs == C.fs
@test C.freqmin == freqmin
@test C.freqmax == freqmax
@test C.corr == Cnew.corr

# test zerophase
A = rand(T,N,Nwin)
C.corr = A
@test bandpass(C,freqmin,freqmax,corners=corners,zerophase=true).corr != bandpass(C,freqmin,
               freqmax,corners=corners,zerophase=false).corr

# test corners
@test bandpass(C,freqmin,freqmax,corners=4,zerophase=true).corr != bandpass(C,freqmin,
              freqmax,corners=3,zerophase=true).corr

# test filtering SeisChannel
C = SeisChannel()
C.x = rand(T,N)
C.fs = fs
C.t = [1 starttime * 1e6;length(C.x) 0]
Cnew = bandpass(C,freqmin,freqmax)
@test Cnew.fs == fs
bandpass!(C,freqmin,freqmax)
@test C.fs == C.fs
@test C.x == Cnew.x

# test zerophase -> there is no option for zerophase here
A = rand(T,N)
C.x = A
@test bandpass(C,freqmin,freqmax,corners=corners,zerophase=true).x == bandpass(C,freqmin,
               freqmax,corners=corners,zerophase=false).x

# test corners
@test bandpass(C,freqmin,freqmax,corners=4,zerophase=true).x != bandpass(C,freqmin,
              freqmax,corners=3,zerophase=true).x

# test filtering SeisData
S = SeisData(2)
S.x[1] = rand(T,N)
S.fs[1] = fs
S.t[1] = [1 starttime * 1e6;length(C.x) 0]
S.x[2] = rand(T,N)
S.fs[2] = fs
S.t[2] = [1 starttime * 1e6;length(C.x) 0]
Snew = bandpass(S,freqmin,freqmax)
@test Snew.fs == S.fs
bandpass!(S,freqmin,freqmax)
@test all(S.fs .== fs)
@test S.x[1] == Snew.x[1]
@test S.x[2] == Snew.x[2]

# test zerophase -> there is no option for zerophase here
A = rand(T,N,Nwin)
S.x[1] = rand(T,N)
S.x[2] = rand(T,N)
@test bandpass(S,freqmin,freqmax,corners=corners,zerophase=true).x[1] == bandpass(S,freqmin,
               freqmax,corners=corners,zerophase=false).x[1]
@test bandpass(S,freqmin,freqmax,corners=corners,zerophase=true).x[2] == bandpass(S,freqmin,
              freqmax,corners=corners,zerophase=false).x[2]

# test corners
@test bandpass(S,freqmin,freqmax,corners=4,zerophase=true).x[1] != bandpass(S,freqmin,
              freqmax,corners=3,zerophase=true).x[1]
@test bandpass(S,freqmin,freqmax,corners=4,zerophase=true).x[2] != bandpass(S,freqmin,
            freqmax,corners=3,zerophase=true).x[2]

# NEEDED: test GPUArray
end

@testset "bandstop" begin

# create white noise
A = rand(T,N,Nwin)
B = bandstop(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
# check size and type stability
@test size(B) == size(A)
@test eltype(B) == eltype(A)
bandstop!(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
@test A == B

# test power in pass bands
freqs = FFTW.rfftfreq(N,fs)
lowband = findall(0 + eps(T) .< freqs .< freqmin)
passband = findall(freqmin .<= freqs .<= freqmax)
highband = findall(freqs .> freqmax)
F = rfft(A,1)
Fpow = abs.(F)
# test lowband
@test all(median(Fpow[lowband,:],dims=1) .> median(Fpow[passband,:],dims=1))

# test highband
@test all(median(Fpow[highband,:],dims=1) .> median(Fpow[passband,:],dims=1))

# test zerophase
A = rand(T,N,Nwin)
@test bandstop(A,freqmin,freqmax,fs,corners=corners,zerophase=true) != bandstop(A,freqmin,
               freqmax,fs,corners=corners,zerophase=false)

# test corners
@test bandstop(A,freqmin,freqmax,fs,corners=4,zerophase=true) != bandstop(A,freqmin,
               freqmax,fs,corners=3,zerophase=true)

# test RawData filtering
R = RawData()
R.x = A
R.fs = fs
Rnew = bandstop(R,freqmin,freqmax)
@test Rnew.fs == R.fs
bandstop!(R,freqmin,freqmax)
@test R.fs == R.fs
@test R.x == Rnew.x

# test zerophase
A = rand(T,N,Nwin)
R.x = A
@test bandstop(R,freqmin,freqmax,corners=corners,zerophase=true).x != bandstop(R,freqmin,
               freqmax,corners=corners,zerophase=false).x

# test corners
@test bandstop(R,freqmin,freqmax,corners=4,zerophase=true).x != bandstop(R,freqmin,
              freqmax,corners=3,zerophase=true).x

# test CorrData filtering
C = CorrData()
C.corr = A
C.fs = fs
Cnew = bandstop(C,freqmin,freqmax)
@test Cnew.fs == C.fs
bandstop!(C,freqmin,freqmax)
@test C.fs == C.fs
@test C.corr == Cnew.corr

# test zerophase
A = rand(T,N,Nwin)
C.corr = A
@test bandstop(C,freqmin,freqmax,corners=corners,zerophase=true).corr != bandstop(C,freqmin,
               freqmax,corners=corners,zerophase=false).corr

# test corners
@test bandstop(C,freqmin,freqmax,corners=4,zerophase=true).corr != bandstop(C,freqmin,
              freqmax,corners=3,zerophase=true).corr

# test filtering SeisChannel
C = SeisChannel()
C.x = rand(T,N)
C.fs = fs
C.t = [1 starttime * 1e6;length(C.x) 0]
Cnew = bandstop(C,freqmin,freqmax)
@test Cnew.fs == fs
bandstop!(C,freqmin,freqmax)
@test C.fs == C.fs
@test C.x == Cnew.x

# test zerophase -> there is no option for zerophase here
A = rand(T,N)
C.x = A
@test bandstop(C,freqmin,freqmax,corners=corners,zerophase=true).x == bandstop(C,freqmin,
               freqmax,corners=corners,zerophase=false).x

# test corners
@test bandstop(C,freqmin,freqmax,corners=4,zerophase=true).x != bandpass(C,freqmin,
              freqmax,corners=3,zerophase=true).x

# test filtering SeisData
S = SeisData(2)
S.x[1] = rand(T,N)
S.fs[1] = fs
S.t[1] = [1 starttime * 1e6;length(C.x) 0]
S.x[2] = rand(T,N)
S.fs[2] = fs
S.t[2] = [1 starttime * 1e6;length(C.x) 0]
Snew = bandstop(S,freqmin,freqmax)
@test Snew.fs == S.fs
bandstop!(S,freqmin,freqmax)
@test all(S.fs .== fs)
@test S.x[1] == Snew.x[1]
@test S.x[2] == Snew.x[2]

# test zerophase -> there is no option for zerophase here
A = rand(T,N,Nwin)
S.x[1] = rand(T,N)
S.x[2] = rand(T,N)
@test bandstop(S,freqmin,freqmax,corners=corners,zerophase=true).x[1] == bandstop(S,freqmin,
               freqmax,corners=corners,zerophase=false).x[1]
@test bandstop(S,freqmin,freqmax,corners=corners,zerophase=true).x[2] == bandstop(S,freqmin,
              freqmax,corners=corners,zerophase=false).x[2]

# test corners
@test bandstop(S,freqmin,freqmax,corners=4,zerophase=true).x[1] != bandstop(S,freqmin,
              freqmax,corners=3,zerophase=true).x[1]
@test bandstop(S,freqmin,freqmax,corners=4,zerophase=true).x[2] != bandstop(S,freqmin,
            freqmax,corners=3,zerophase=true).x[2]

# NEEDED: test GPUArray
end

# test bandpass
@testset "lowpass" begin

# create white noise
A = rand(T,N,Nwin)
B = lowpass(A,freqmax,fs,corners=corners,zerophase=zerophase)
# check size and type stability
@test size(B) == size(A)
@test eltype(B) == eltype(A)
lowpass!(A,freqmax,fs,corners=corners,zerophase=zerophase)
@test A == B

# test power in pass bands
freqs = FFTW.rfftfreq(N,fs)
lowband = findall(0 + eps(T) .< freqs .< freqmax)
highband = findall(freqs .>= freqmax)
F = rfft(A,1)
Fpow = abs.(F)
# test filter
@test all(median(Fpow[lowband,:],dims=1) .> median(Fpow[highband,:],dims=1))

# test zerophase
A = rand(T,N,Nwin)
@test lowpass(A,freqmax,fs,corners=corners,zerophase=true) != lowpass(A,
               freqmax,fs,corners=corners,zerophase=false)

# test corners
@test lowpass(A,freqmax,fs,corners=4,zerophase=true) != lowpass(A,
               freqmax,fs,corners=3,zerophase=true)

# test RawData filtering
R = RawData()
R.x = A
R.fs = fs
Rnew = lowpass(R,freqmax)
@test Rnew.fs == R.fs
@test Rnew.freqmax == freqmax
lowpass!(R,freqmax)
@test R.fs == R.fs
@test R.freqmax == freqmax
@test R.x == Rnew.x

# test zerophase
A = rand(T,N,Nwin)
R.x = A
@test lowpass(R,freqmax,corners=corners,zerophase=true).x != lowpass(R,
               freqmax,corners=corners,zerophase=false).x

# test corners
@test lowpass(R,freqmax,corners=4,zerophase=true).x != lowpass(R,
              freqmax,corners=3,zerophase=true).x

# test CorrData filtering
C = CorrData()
C.corr = A
C.fs = fs
Cnew = lowpass(C,freqmax)
@test Cnew.fs == C.fs
@test Cnew.freqmax == freqmax
lowpass!(C,freqmax)
@test C.fs == C.fs
@test C.freqmax == freqmax
@test C.corr == Cnew.corr

# test zerophase
A = rand(T,N,Nwin)
C.corr = A
@test lowpass(C,freqmax,corners=corners,zerophase=true).corr != lowpass(C,
               freqmax,corners=corners,zerophase=false).corr

# test corners
@test lowpass(C,freqmax,corners=4,zerophase=true).corr != lowpass(C,
              freqmax,corners=3,zerophase=true).corr

# test filtering SeisChannel
C = SeisChannel()
C.x = rand(T,N)
C.fs = fs
C.t = [1 starttime * 1e6;length(C.x) 0]
Cnew = lowpass(C,freqmax)
@test Cnew.fs == fs
lowpass!(C,freqmax)
@test C.fs == C.fs
@test C.x == Cnew.x

# test zerophase -> there is no option for zerophase here
A = rand(T,N)
C.x = A
@test lowpass(C,freqmax,corners=corners,zerophase=true).x == lowpass(C,
               freqmax,corners=corners,zerophase=false).x

# test corners
@test lowpass(C,freqmax,corners=4,zerophase=true).x != lowpass(C,
              freqmax,corners=3,zerophase=true).x

# test filtering SeisData
S = SeisData(2)
S.x[1] = rand(T,N)
S.fs[1] = fs
S.t[1] = [1 starttime * 1e6;length(C.x) 0]
S.x[2] = rand(T,N)
S.fs[2] = fs
S.t[2] = [1 starttime * 1e6;length(C.x) 0]
Snew = lowpass(S,freqmax)
@test Snew.fs == S.fs
lowpass!(S,freqmax)
@test all(S.fs .== fs)
@test S.x[1] == Snew.x[1]
@test S.x[2] == Snew.x[2]

# test zerophase -> there is no option for zerophase here
A = rand(T,N,Nwin)
S.x[1] = rand(T,N)
S.x[2] = rand(T,N)
@test lowpass(S,freqmax,corners=corners,zerophase=true).x[1] == lowpass(S,
               freqmax,corners=corners,zerophase=false).x[1]
@test lowpass(S,freqmax,corners=corners,zerophase=true).x[2] == lowpass(S,
              freqmax,corners=corners,zerophase=false).x[2]

# test corners
@test lowpass(S,freqmax,corners=4,zerophase=true).x[1] != lowpass(S,
              freqmax,corners=3,zerophase=true).x[1]
@test lowpass(S,freqmax,corners=4,zerophase=true).x[2] != lowpass(S,
            freqmax,corners=3,zerophase=true).x[2]

# NEEDED: test GPUArray
end


@testset "highpass" begin

# create white noise
A = rand(T,N,Nwin)
B = highpass(A,freqmin,fs,corners=corners,zerophase=zerophase)
# check size and type stability
@test size(B) == size(A)
@test eltype(B) == eltype(A)
highpass!(A,freqmin,fs,corners=corners,zerophase=zerophase)
@test A == B

# test power in pass bands
freqs = FFTW.rfftfreq(N,fs)
lowband = findall(0 + eps(T) .< freqs .< freqmax)
highband = findall(freqs .>= freqmax)
F = rfft(A,1)
Fpow = abs.(F)
# test filter
@test all(median(Fpow[lowband,:],dims=1) .< median(Fpow[highband,:],dims=1))

# test zerophase
A = rand(T,N,Nwin)
@test highpass(A,freqmin,fs,corners=corners,zerophase=true) != highpass(A,
               freqmin,fs,corners=corners,zerophase=false)

# test corners
@test highpass(A,freqmin,fs,corners=4,zerophase=true) != highpass(A,
               freqmin,fs,corners=3,zerophase=true)

# test RawData filtering
R = RawData()
R.x = A
R.fs = fs
Rnew = highpass(R,freqmin)
@test Rnew.fs == R.fs
@test Rnew.freqmin == freqmin
highpass!(R,freqmin)
@test R.fs == R.fs
@test R.freqmin == freqmin
@test R.x == Rnew.x

# test zerophase
A = rand(T,N,Nwin)
R.x = A
@test highpass(R,freqmin,corners=corners,zerophase=true).x != highpass(R,
               freqmin,corners=corners,zerophase=false).x

# test corners
@test highpass(R,freqmin,corners=4,zerophase=true).x != highpass(R,
              freqmin,corners=3,zerophase=true).x

# test CorrData filtering
C = CorrData()
C.corr = A
C.fs = fs
Cnew = highpass(C,freqmin)
@test Cnew.fs == C.fs
@test Cnew.freqmin == freqmin
highpass!(C,freqmin)
@test C.fs == C.fs
@test C.freqmin == freqmin
@test C.corr == Cnew.corr

# test zerophase
A = rand(T,N,Nwin)
C.corr = A
@test highpass(C,freqmin,corners=corners,zerophase=true).corr != highpass(C,
               freqmin,corners=corners,zerophase=false).corr

# test corners
@test highpass(C,freqmin,corners=4,zerophase=true).corr != highpass(C,
              freqmin,corners=3,zerophase=true).corr

# test filtering SeisChannel
C = SeisChannel()
C.x = rand(T,N)
C.fs = fs
C.t = [1 starttime * 1e6;length(C.x) 0]
Cnew = highpass(C,freqmin)
@test Cnew.fs == fs
highpass!(C,freqmin)
@test C.fs == C.fs
@test C.x == Cnew.x

# test zerophase -> there is no option for zerophase here
A = rand(T,N)
C.x = A
@test highpass(C,freqmin,corners=corners,zerophase=true).x == highpass(C,
               freqmin,corners=corners,zerophase=false).x

# test corners
@test highpass(C,freqmin,corners=4,zerophase=true).x != highpass(C,
              freqmin,corners=3,zerophase=true).x

# test filtering SeisData
S = SeisData(2)
S.x[1] = rand(T,N)
S.fs[1] = fs
S.t[1] = [1 starttime * 1e6;length(C.x) 0]
S.x[2] = rand(T,N)
S.fs[2] = fs
S.t[2] = [1 starttime * 1e6;length(C.x) 0]
Snew = highpass(S,freqmin)
@test Snew.fs == S.fs
highpass!(S,freqmin)
@test all(S.fs .== fs)
@test S.x[1] == Snew.x[1]
@test S.x[2] == Snew.x[2]

# test zerophase -> there is no option for zerophase here
A = rand(T,N,Nwin)
S.x[1] = rand(T,N)
S.x[2] = rand(T,N)
@test highpass(S,freqmin,corners=corners,zerophase=true).x[1] == highpass(S,
               freqmin,corners=corners,zerophase=false).x[1]
@test highpass(S,freqmin,corners=corners,zerophase=true).x[2] == highpass(S,
              freqmin,corners=corners,zerophase=false).x[2]

# test corners
@test highpass(S,freqmin,corners=4,zerophase=true).x[1] != highpass(S,
              freqmin,corners=3,zerophase=true).x[1]
@test highpass(S,freqmin,corners=4,zerophase=true).x[2] != highpass(S,
            freqmin,corners=3,zerophase=true).x[2]

# NEEDED: test GPUArray
end

# test envelope
@testset "envelope" begin
t = collect(0.:π/100:10π)
t = T.(t)
A = sin.(t) .* cos.(5 .* t)
up, down = envelope(A)
@test size(up) == size(A)
@test size(down) == size(A)
@test eltype(up) == eltype(A)
@test eltype(down) == eltype(A)
@test all(up .- A .> -1e-6)
@test all(A .- down .> -1e-6)

end
