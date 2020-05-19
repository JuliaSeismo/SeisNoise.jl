# test phase shifting

T = Float32
fs = 100.  # sampling rate
N = 2^14   # number of points
starttime = d2u(DateTime(Date(now())))

@testset "phase shift" begin

C = SeisChannel()
C.x = rand(T,N)
C.fs = fs
C.t = [1 starttime * 1e6;length(C.x) 0]

# test no shift
Cnew = phase_shift(C,ϕshift=true)
@test Cnew == C

# test phase shift + time shift
C.t = [1 (starttime + 0.007) * 1e6;length(C.x) 0]
Cnew = phase_shift(C,ϕshift=true)
@test Cnew != C
@test Cnew.t != C.t
@test Cnew.x != C.x

# test just time shift
Cnew = phase_shift(C,ϕshift=false)
@test Cnew != C
@test Cnew.t != C.t
@test Cnew.x == C.x

# test in-place
Cnew = phase_shift(C,ϕshift=true)
phase_shift!(C)
@test Cnew == C
@test Cnew.t == C.t
@test Cnew.x == C.x

## test SeisData
S = SeisData(2)
S.x[1] = rand(T,N)
S.fs[1] = fs
# unchanged starttime
S.t[1] = [1 starttime * 1e6;length(C.x) 0]
S.x[2] = rand(T,N)
S.fs[2] = fs
# offset starttime
S.t[2] = [1 (starttime + 0.007) * 1e6;length(C.x) 0]

# test phase shift + time shift
Snew = phase_shift(S,ϕshift=true)
@test Snew != S
@test Snew[1] == S[1]
@test Snew[2] != S[2]

# test just time shift
Snew = phase_shift(S,ϕshift=false)
@test Snew != S
@test Snew.t[1] == S.t[1]
@test Snew.x[1] == S.x[1]
@test Snew.t[2] != S.t[2]
@test Snew.x[2] == S.x[2]

# test in-place
Snew = phase_shift(S,ϕshift=true)
phase_shift!(S)
@test Snew == S
@test Snew.t[1] == S.t[1]
@test Snew.x[1] == S.x[1]
@test Snew.t[2] == S.t[2]
@test Snew.x[2] == S.x[2]
end
