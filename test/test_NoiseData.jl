# test creation of RawData, FFTData, CorrData

# create RawData from SeisData
fs = 100.
nx = Int(86400 * fs)
S = SeisData(2)
S.id[1] = "N1.SS1.00.HHZ"
S.name[1] = "New Station"
S.loc[1] = GeoLoc()
S.fs[1] = 100.0
S.gain[1] = 2.0
S.resp[1] = PZResp([0.0+0.0*im 1.0+0.767*im; 0.0+0.0*im 1.0-0.767*im])
S.units[1] = ""
S.src[1] = ""
S.notes[1] = ["asdfasdf","dfhdfgh"]
S.misc[1] = Dict{String,Any}("P" => 2.0)
S.t[1] = [1 round(Int, d2u(DateTime(Date(now())))*1e6); nx 0]
S.x[1] = rand(Float32,nx)

S.id[2] = "N2.SS2.00.HHZ"
S.name[2] = "Old Station"
S.loc[2] = GeoLoc()
S.fs[2] = fs
S.gain[2] = 22.0
S.resp[2] = PZResp64(z = [0.0+0.0*im, 0.0+0.0*im], p = [1.0+1.0*im, 1.0-1.0*im])
S.units[2] = "ms/2"
S.src[2] = "file"
S.notes[2] = ["0913840183","klnelgng"]
S.misc[2] = Dict{String,Any}("S" => 6.5)
S.t[2] = [1 round(Int, d2u(DateTime(Date(now())))*1e6); nx 0]
S.x[2] = rand(Float32,nx)

# test empty allocator
@test !isa(try RawData() catch ex ex end, Exception)    # Passes
R = RawData()
@test isempty(R)

# check SeisData to RawData with similar cc_len, cc_step
cc_len = cc_step = 1800.
Ch = S[1]
R = RawData(Ch,cc_len,cc_step)
@test !isa(try RawData(S,cc_len,cc_step) catch ex ex end, Exception)    # Passes
@test length(R.t) == size(R.x,2)        # check number of windows
@test size(R.x,1) == R.cc_len * R.fs    # check length of window
@test R.freqmin == 1 / R.cc_len         # check freqmin
@test R.freqmax == R.fs / 2             # check freqmax
@test ndims(R.x) == 2                   # check dimensions
@test size(R.x,2) == 48                 # check number of windows
@test R == R                            # test equality
@test u2d(R.t[1]) == u2d(Ch.t[1,2] * 1e-6) # test starttimes equal
@test eltype(R.x) == eltype(Ch.x)

# check uncommon cc_len & cc_step
cc_len = π
@test_throws DomainError RawData(Ch,cc_len,cc_step)
cc_len = 1800.
cc_step = π
@test_throws DomainError RawData(Ch,cc_len,cc_step)

# check huge cc_len
cc_len = 86401.
cc_step = 1800.
@test_throws DomainError RawData(Ch,cc_len,cc_step)

# check SeisChannel to RawData with different cc_len, cc_step
cc_len = 1800.
cc_step = 450.
Ch = S[1]
R = RawData(Ch,cc_len,cc_step)
@test !isa(try RawData(Ch,cc_len,cc_step) catch ex ex end, Exception)    # Passes
@test length(R.t) == size(R.x,2)        # check number of windows
@test size(R.x,1) == R.cc_len * R.fs    # check length of window
@test R.freqmin ==  1 / R.cc_len        # check freqmin
@test R.freqmax == R.fs / 2             # check freqmax
@test ndims(R.x) == 2                   # check dimensions
@test size(R.x,2) == 189                # check number of windows
@test R == R                            # test equality
@test u2d(R.t[1]) == u2d(Ch.t[1,2] * 1e-6) # test starttimes equal
@test eltype(R.x) == eltype(Ch.x)

# check uncommon cc_len & cc_step
cc_len = π
@test_throws DomainError RawData(Ch,cc_len,cc_step)
cc_len = 1800.
cc_step = π
@test_throws DomainError RawData(Ch,cc_len,cc_step)

# check huge cc_len
cc_len = 86401.
cc_step = 1800.
@test_throws DomainError RawData(Ch,cc_len,cc_step)

# check SeisChannel to RawData with similar cc_len, cc_step
cc_len = cc_step = 1800.
R = RawData(Ch,cc_len,cc_step)
@test !isa(try RawData(Ch,cc_len,cc_step) catch ex ex end, Exception)    # Passes
@test length(R.t) == size(R.x,2)        # check number of windows
@test size(R.x,1) == R.cc_len * R.fs    # check length of window
@test R.freqmin == 1 / R.cc_len         # check freqmin
@test R.freqmax == R.fs / 2             # check freqmax
@test ndims(R.x) == 2                   # check dimensions
@test size(R.x,2) == 48                 # check number of windows
@test R == R                            # test equality
@test u2d(R.t[1]) == u2d(S.t[1][1,2] * 1e-6) # test starttimes equal
@test eltype(R.x) == eltype(S.x[1])

# check uncommon cc_len & cc_step
cc_len = π
@test_throws DomainError RawData(S,cc_len,cc_step)
cc_len = 1800.
cc_step = π
@test_throws DomainError RawData(S,cc_len,cc_step)

# check huge cc_len
cc_len = 86401.
cc_step = 1800.
@test_throws DomainError RawData(S,cc_len,cc_step)

# check SeisData to RawData with different cc_len, cc_step
cc_len = 1800.
cc_step = 450.
R = RawData(S,cc_len,cc_step)
@test !isa(try RawData(S,cc_len,cc_step) catch ex ex end, Exception)    # Passes
@test length(R.t) == size(R.x,2)        # check number of windows
@test size(R.x,1) == R.cc_len * R.fs    # check length of window
@test R.freqmin ==  1 / R.cc_len        # check freqmin
@test R.freqmax == R.fs / 2             # check freqmax
@test size(R.x,2) == 189                # check number of windows
@test ndims(R.x) == 2                   # check dimensions
@test R == R                            # test equality
@test u2d(R.t[1]) == u2d(S.t[1][1,2] * 1e-6) # test starttimes equal
@test eltype(R.x) == eltype(S.x[1])

# check uncommon cc_len & cc_step
cc_len = π
@test_throws DomainError RawData(S,cc_len,cc_step)
cc_len = 1800.
cc_step = π
@test_throws DomainError RawData(S,cc_len,cc_step)

# check huge cc_len
cc_len = 86401.
cc_step = 1800.
@test_throws DomainError RawData(S,cc_len,cc_step)

## test when starttime is not aligned with cc_len
# this SeisData starts at 00:16:40
S.t[1] = [1 round(Int, (d2u(DateTime(Date(now()))) + 1e3) *1e6 ); nx 0]
cc_len = cc_step = 1800.
Ch = S[1]
R = RawData(S,cc_len,cc_step)
@test !isa(try RawData(Ch,cc_len,cc_step) catch ex ex end, Exception)    # Passes
@test length(R.t) == size(R.x,2)        # check number of windows
@test size(R.x,1) == R.cc_len * R.fs    # check length of window
@test R.freqmin == 1 / R.cc_len         # check freqmin
@test R.freqmax == R.fs / 2             # check freqmax
@test ndims(R.x) == 2                   # check dimensions
@test size(R.x,2) == 47                 # check number of windows
@test eltype(R.x) == eltype(S.x[1])

# test CPU/ GPU
@test R == cpu(R)
@test cpu(R) == cpu(R)
@test R == cpu(gpu(R))

# test in
@test in(R.id,R)

# test addition
R = RawData(S[1],cc_len,cc_step)
R1 = RawData()
@test append!(R1,R) == R
@test append!(R1,R) == R + R

## test FFTData

# test empty allocator
@test !isa(try FFTData() catch ex ex end, Exception)    # Passes
F = FFTData()
@test isempty(F)

# test FFTData creation
diveven(x) = Int(ceil(x / 2)) + iseven(Int(x)) # function to check rfft size
S.t[1] = [1 round(Int, d2u(DateTime(Date(now())))*1e6); nx 0]
cc_len = cc_step = 1800.
R = RawData(S,cc_len,cc_step)
FFT = rfft(R.x,1)

@test !isa(try FFTData(R.name, R.id, R.loc, R.fs, R.gain, R.freqmin,
                       R.freqmax, R.cc_len, R.cc_step, false, "", R.resp,
                       R.misc, R.notes, R.t, FFT) catch ex ex end, Exception)    # Passes
F = FFTData(R.name, R.id, R.loc, R.fs, R.gain, R.freqmin,
                       R.freqmax, R.cc_len, R.cc_step, false, "", R.resp,
                       R.misc, R.notes, R.t, FFT)
@test length(F.t) == size(F.fft,2)              # check number of windows
@test size(F.fft,1) == diveven(F.cc_len * F.fs) # check length of window
@test F.freqmin ==  1 / F.cc_len                # check freqmin
@test F.freqmax == F.fs / 2                     # check freqmax
@test size(F.fft,2) == 48                       # check number of windows
@test ndims(F.fft) == 2                         # check dimensions
@test F == F                                    # test equality
@test u2d(F.t[1]) == u2d(S.t[1][1,2] * 1e-6)    # test starttimes equal
@test eltype(F.fft) <: Complex

# test CPU/ GPU
@test F == cpu(F)
@test cpu(F) == cpu(F)
@test F == cpu(gpu(F))

# test in
@test in(F.id,F)

# test addition
F1 = FFTData()
@test append!(F1,F) == F
@test append!(F1,F) == F + F

## test CorrData
# test empty allocator
@test !isa(try CorrData() catch ex ex end, Exception)    # Passes
C = CorrData()
@test isempty(C)

# test CorrData creation
maxlag = 100.
cc_len = cc_step = 1800.
R1 = RawData(S[1],cc_len,cc_step)
R2 = RawData(S[2],cc_len,cc_step)
F1 = rfft(R1)
F2 = rfft(R2)
comp = F1.name[end] * F2.name[end]
corr = correlate(F1.fft,F2.fft,Int(F1.cc_len * F1.fs),Int(maxlag*F1.fs))

@test !isa(try CorrData(F1.name, F1.id, F1.loc, comp, false, "CC", F1.fs, F1.gain,
  F1.freqmin, F1.freqmax, F1.cc_len, F1.cc_step, F1.whitened, F1.time_norm,
  F1.resp, F1.misc, F1.notes, 0.,0., 0., maxlag, F1.t, corr) catch ex ex end, Exception)    # Passes
C = CorrData(F1.name, F1.id, F1.loc, comp, false, "CC", F1.fs, F1.gain,
  F1.freqmin, F1.freqmax, F1.cc_len, F1.cc_step, F1.whitened, F1.time_norm,
  F1.resp, F1.misc, F1.notes, 0.,0., 0., maxlag, F1.t, corr)
@test length(C.t) == size(C.corr,2)             # check number of windows
@test size(C.corr,1) == Int(2 * C.maxlag * C.fs + 1) # check length of window
@test C.freqmin ==  1 / C.cc_len                # check freqmin
@test C.freqmax == C.fs / 2                     # check freqmax
@test size(C.corr,2) == 48                      # check number of windows
@test ndims(C.corr) == 2                        # check dimensions
@test C == C                                    # test equality
@test u2d(C.t[1]) == u2d(S.t[1][1,2] * 1e-6)    # test starttimes equal
@test eltype(C.corr) <: Real

# test CPU/ GPU
@test C == cpu(C)
@test cpu(C) == cpu(C)
@test C == cpu(gpu(C))

# test in
@test in(C.id,C)

# test addition
C1 = CorrData()
@test append!(C1,C) == C
@test append!(C1,C) == C + C
