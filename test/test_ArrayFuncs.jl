# test array functions
N = 10000 # data length

# function to test similarity of vectors
cosinedist(a,b) = sum(a .* b, dims=1) ./ sqrt.(sum(a .^ 2, dims=1)) ./ sqrt.(sum(b .^ 2,dims=1))

## 1D detrending
@testset "Test detrending" begin
A = rand(N)
A .-= SeisNoise.mean(A) # remove mean
x = collect(1:N)
trendA = A .+ x         # add trend

# test 1D array in-place
@test !isa(try detrend!(trendA) catch ex ex end, Exception)    # Passes
@test cosinedist(A,trendA)[1]  > 0.999
@test eltype(trendA) == eltype(A) # check type stability
@test size(trendA) == size(A)     # check size remains same

# test allocating 1D
trendB = detrend(A .+ x)
@test cosinedist(A,trendB)[1]  > 0.999
@test eltype(trendB) == eltype(A) # check type stability
@test size(trendB) == size(A)     # check size remains same

## 2D detrending
A = rand(N,10) # detrend 10 columns
A .-= SeisNoise.mean(A,dims=1) # remove mean
trendA = A .+ x

# test 2D array in-place
@test !isa(try detrend!(trendA) catch ex ex end, Exception)    # Passes
@test all(cosinedist(A,trendA)  .> 0.999)

# test allocating 2D
trendB = detrend(A .+ x)
@test all(cosinedist(A,trendB)  .> 0.999)

## test RawData
# test in-place
R = RawData()
A = rand(N,10) # detrend 10 columns
A .-= SeisNoise.mean(A,dims=1) # remove mean
R.x = A .+ x
@test !isa(try detrend!(R) catch ex ex end, Exception)    # Passes
@test isa(detrend!(R),Nothing)
@test all(cosinedist(A,R.x)  .> 0.999)

# test allocating
R = RawData()
A = rand(N,10) # detrend 10 columns
A .-= SeisNoise.mean(A,dims=1) # remove mean
R.x = A .+ x
@test !isa(try detrend(R) catch ex ex end, Exception)    # Passes
@test isa(detrend(R),RawData)
Rnew = detrend(R)
@test all(cosinedist(A,Rnew.x)  .> 0.999)

## test CorrData
# test in-place
C = CorrData()
A = rand(N,10) # detrend 10 columns
A .-= SeisNoise.mean(A,dims=1) # remove mean
C.corr = A .+ x
@test !isa(try detrend!(C) catch ex ex end, Exception)    # Passes
@test isa(detrend!(C),Nothing)
@test all(cosinedist(A,C.corr)  .> 0.999)

# test allocating
C = CorrData()
A = rand(N,10) # detrend 10 columns
A .-= SeisNoise.mean(A,dims=1) # remove mean
C.corr = A .+ x
@test !isa(try detrend(C) catch ex ex end, Exception)    # Passes
@test isa(detrend(C),CorrData)
Cnew = detrend(C)
@test all(cosinedist(A,Cnew.corr)  .> 0.999)

## test failure on Integer data
A = rand(-2^30 : 2^30, N)
@test_throws MethodError detrend(A)
end

## test demean
@testset "Test demeaning" begin
# 1D demean in-place
A = Float64.(collect(-N:N))
x = rand() * N
meanA = A .+ x
@test !isa(try demean!(meanA) catch ex ex end, Exception)    # Passes
@test isapprox(A,meanA) # test all close
@test eltype(meanA) == eltype(A)
@test size(meanA) == size(A)

# test allocating 1D
meanB = demean(A .+ x)
@test isapprox(A,meanB) # test all close

# 2D demean
A = repeat(collect(-N:N),outer=(1,10))
x = rand(1,10) * N
meanA = A .+ x
@test !isa(try demean!(meanA) catch ex ex end, Exception)    # Passes
@test isapprox(A,meanA) # test all close

# test allocating 2D
meanB = demean(A .+ x)
@test isapprox(A,meanB) # test all close

# test failure on integer data
A = rand(-2^30 : 2^30, N)
@test_throws MethodError demean(A)

## test RawData
# test in-place
R = RawData()
A = repeat(collect(-N:N),outer=(1,10))
x = rand(1,10) * N
R.x = A .+ x
@test !isa(try demean!(R) catch ex ex end, Exception)    # Passes
@test isa(demean!(R),Nothing)
@test isapprox(A,R.x)

# test allocating
R = RawData()
R.x = A .+ x
@test !isa(try detrend(R) catch ex ex end, Exception)    # Passes
@test isa(detrend(R),RawData)
Rnew = demean(R)
@test isapprox(A,Rnew.x)

## test CorrData
# test in-place
C = CorrData()
C.corr = A .+ x
@test !isa(try demean!(C) catch ex ex end, Exception)    # Passes
@test isa(demean!(C),Nothing)
@test isapprox(A,C.corr)

# test allocating
C = CorrData()
C.corr = A .+ x
@test !isa(try demean(C) catch ex ex end, Exception)    # Passes
@test isa(detrend(C),CorrData)
Cnew = demean(C)
@test isapprox(A,Cnew.corr)
end

## test taper
@testset "Test tapering" begin
# test 1D in-place
A = ones(N)
taperA = deepcopy(A)
max_percentage = 0.05
fs = 100.
max_length = 20.

# check default taper lengths
Nperc = Int(max_percentage * N)
Nlength = Int(max_length * fs)

taperA = deepcopy(A)
@test !isa(try taper!(taperA,fs) catch ex ex end, Exception)
@test all(taperA[1:min(Nperc,Nlength)] .!= 1)
@test all(taperA[end-min(Nperc,Nlength):end] .!= 1)
@test all(taperA[min(Nperc,Nlength)+1:end-min(Nperc,Nlength)-1] .== 1)
@test eltype(taperA) == eltype(A) # check type stability
@test size(taperA) == size(A)     # check size remains same

# test 1D allocating
@test !isa(try taper(A,fs) catch ex ex end, Exception)
@test isa(taper(A,fs),typeof(A))
taperB = taper(A,fs)
@test all(taperA .== taperB)
@test eltype(taperB) == eltype(A) # check type stability

# change max_percentage
max_percentage = 0.01
Nperc = Int(max_percentage * N)
taperA = deepcopy(A)
taper!(taperA,fs,max_percentage=max_percentage)
@test all(taperA[1:Nperc] .!= 1)
@test all(taperA[end-Nperc:end] .!= 1)
@test all(taperA[Nperc+1:end-Nperc-1] .== 1)

# change max_percentage
max_length = 1.
Nlength = Int(max_length * fs)
taperA = deepcopy(A)
taper!(taperA,fs,max_length=max_length)
@test all(taperA[1:Nlength] .!= 1)
@test all(taperA[end-Nlength:end] .!= 1)
@test all(taperA[Nlength+1:end-Nlength-1] .== 1)

# test 2D in-place
A = ones(N,10)
taperA = deepcopy(A)
max_percentage = 0.05
fs = 100.
max_length = 20.

# check default taper lengths
Nperc = Int(max_percentage * N)
Nlength = Int(max_length * fs)

taperA = deepcopy(A)
@test !isa(try taper!(taperA,fs) catch ex ex end, Exception)
@test all(taperA[1:min(Nperc,Nlength),:] .!= 1)
@test all(taperA[end-min(Nperc,Nlength):end,:] .!= 1)
@test all(taperA[min(Nperc,Nlength)+1:end-min(Nperc,Nlength)-1,:] .== 1)

# change max_percentage
max_percentage = 0.01
Nperc = Int(max_percentage * N)
taperA = deepcopy(A)
taper!(taperA,fs,max_percentage=max_percentage)
@test all(taperA[1:Nperc,:] .!= 1)
@test all(taperA[end-Nperc:end,:] .!= 1)
@test all(taperA[Nperc+1:end-Nperc-1,:] .== 1)

# change max_percentage
max_length = 1.
Nlength = Int(max_length * fs)
taperA = deepcopy(A)
taper!(taperA,fs,max_length=max_length)
@test all(taperA[1:Nlength,:] .!= 1)
@test all(taperA[end-Nlength:end,:] .!= 1)
@test all(taperA[Nlength+1:end-Nlength-1,:] .== 1)


# test using Real for fs and max_length
taperA = deepcopy(A)
taper!(taperA,Int(fs),max_length=Int(max_length))
@test all(taperA[1:Nlength,:] .!= 1)
@test all(taperA[end-Nlength:end,:] .!= 1)
@test all(taperA[Nlength+1:end-Nlength-1,:] .== 1)

# test RawData
# test in-place
R = RawData()
A = ones(N,10) # detrend 10 columns
max_percentage = 0.05
fs = 100.
max_length = 20.
Nperc = Int(max_percentage * N)
Nlength = Int(max_length * fs)
R.x = A
R.fs = fs
@test !isa(try taper!(R) catch ex ex end, Exception)    # Passes
@test isa(taper!(R),Nothing)
@test all(R.x[1:min(Nperc,Nlength),:] .!= 1)
@test all(R.x[end-min(Nperc,Nlength):end,:] .!= 1)
@test all(R.x[min(Nperc,Nlength)+1:end-min(Nperc,Nlength)-1,:] .== 1)

# test allocating
R = RawData()
R.fs = fs
R.x = A
@test !isa(try taper(R) catch ex ex end, Exception)    # Passes
@test isa(taper(R),RawData)
Rnew = taper(R)
@test all(Rnew.x[1:min(Nperc,Nlength),:] .!= 1)
@test all(Rnew.x[end-min(Nperc,Nlength):end,:] .!= 1)
@test all(Rnew.x[min(Nperc,Nlength)+1:end-min(Nperc,Nlength)-1,:] .== 1)
@test all(Rnew.x[1:min(Nperc,Nlength),:] .!= 1)
@test all(Rnew.x[end-min(Nperc,Nlength):end,:] .!= 1)
@test all(Rnew.x[min(Nperc,Nlength)+1:end-min(Nperc,Nlength)-1,:] .== 1)

# test with integers 
Rnew = taper(R,max_length=Int(max_length))

## test CorrData
# test in-place
C = CorrData()
C.fs = fs
C.corr = A
@test !isa(try taper!(C) catch ex ex end, Exception)    # Passes
@test isa(taper!(C),Nothing)
@test all(C.corr[1:min(Nperc,Nlength),:] .!= 1)
@test all(C.corr[end-min(Nperc,Nlength):end,:] .!= 1)
@test all(C.corr[min(Nperc,Nlength)+1:end-min(Nperc,Nlength)-1,:] .== 1)

# test allocating
C = CorrData()
C.fs = fs
C.corr = A
@test !isa(try taper(C) catch ex ex end, Exception)    # Passes
@test isa(taper(C),CorrData)
Cnew = taper(C)
@test all(Cnew.corr[1:min(Nperc,Nlength),:] .!= 1)
@test all(Cnew.corr[end-min(Nperc,Nlength):end,:] .!= 1)
@test all(Cnew.corr[min(Nperc,Nlength)+1:end-min(Nperc,Nlength)-1,:] .== 1)

# test with Integers 
Cnew = taper(C,max_length=Int(max_length))
@test all(Cnew.corr[1:min(Nperc,Nlength),:] .!= 1)
@test all(Cnew.corr[end-min(Nperc,Nlength):end,:] .!= 1)
@test all(Cnew.corr[min(Nperc,Nlength)+1:end-min(Nperc,Nlength)-1,:] .== 1)

## test failure on Integer data
A = ones(Int,N)
@test_throws MethodError taper(A,fs)

## test hanningwindow
A = ones(N)
@test !isa(try hanningwindow(A,2*Nlength) catch ex ex end, Exception)
win = hanningwindow(A,2 * Nlength)
@test eltype(win) == eltype(A)     # test type stability
@test length(win) == 2 * Nlength   # test vector length
@test isapprox(win[1:Nlength],win[end:-1:Nlength+1]) # test taper mirroring
end
