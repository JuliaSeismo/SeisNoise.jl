# test tools for matrix smoothing

T = Float32
@testset "SNR" begin
A = rand(100,100) .+ 10. .* sin.(0.:2π/100:2π) .- 0.5
A .= T.(A)
# test high SNR
@test all(snr(A) .> 10)
# test type stability
@test eltype(snr(A)) == eltype(A)

C = CorrData()
C.corr = A
@test all(snr(C) .> 10)
end

@testset "smoothing" begin
A = rand(Float32,1000,10)
presmooth = mean(abs.(diff(A,dims=1)),dims=1)
el = eltype(A)
# test return type
@test smooth!(A) == nothing
# test type stability
@test eltype(A) == el
# test smoothness
@test all(presmooth .> mean(abs.(diff(A,dims=1)),dims=1))

@test eltype(smooth(A)) == eltype(A)
@test size(smooth(A)) == size(A)
# test smooth with window lendth of 3
B = smooth(A)
@test isapprox(mean(A[3:9,1]),B[6,1])
end

@testset "abs_max" begin
A = (rand(T,1000,10) .- T(0.5)) .* 1000
@test abs_max!(A) == nothing
@test maximum(A) <= 1.
@test minimum(A) >= -1.
A[:,1] .= 0
@test_throws DomainError abs_max!(A)

# test allocating
A = (rand(T,1000,10) .- T(0.5)) .* 1000
@test size(abs_max(A)) == size(A)
@test eltype(abs_max(A)) == eltype(A)
B = abs_max(A)
@test maximum(B) <= 1.
@test minimum(B) >= -1.

# test CorrData
C = CorrData()
A = (rand(T,1000,10) .- T(0.5)) .* 1000
C.corr = A
@test abs_max!(C) == nothing
@test maximum(C.corr) <= 1.
@test minimum(C.corr) >= -1

# test CorrData
C = CorrData()
A = (rand(T,1000,10) .- T(0.5)) .* 1000
C.corr = A
@test isa(abs_max(C),CorrData)
c = abs_max(C)
@test maximum(c.corr) <= 1.
@test minimum(c.corr) >= -1
end

@testset "standardize" begin
A = (rand(T,1000,10) .- T(0.5)) .* 1000
@test standardize!(A) == nothing
@test all(std(A,dims=1) .≈ T(1.))

# test allocating
A = (rand(T,1000,10) .- T(0.5)) .* 1000
@test size(standardize(A)) == size(A)
@test eltype(standardize(A)) == eltype(A)
B = standardize(A)
@test all(std(B,dims=1) .≈ T(1.))

# test CorrData
C = CorrData()
A = (rand(T,1000,10) .- T(0.5)) .* 1000
C.corr = A
@test standardize!(C) == nothing
@test all(std(C.corr,dims=1) .≈ T(1.))

# test CorrData
C = CorrData()
A = (rand(T,1000,10) .- T(0.5)) .* 1000
C.corr = A
@test isa(standardize(C),CorrData)
c = standardize(C)
@test all(std(c.corr,dims=1) .≈ T(1.))
end

@testset "mad" begin
A = rand(100,10)
A .= T.(A)
# test high mad
@test all(mad(A) .< 1)
# test type stability
@test eltype(mad(A)) == eltype(A)
@test size(mad(A),2) .== size(A,2)

C = CorrData()
C.corr = A
@test all(mad(C) .< 1)
end

@testset "std_threshold" begin
A = rand(T,1000,10)
A[100,1] = T(1000)
@test std_threshold(A,5.) == collect(2:10)
end
