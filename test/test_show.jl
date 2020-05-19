# test show
@testset "show" begin
(path, io) = mktemp()
# show
redirect_stdout(io) do
@test !isa(try show(RawData()) catch ex ex end, Exception)
@test !isa(try show(FFTData()) catch ex ex end, Exception)
@test !isa(try show(CorrData()) catch ex ex end, Exception)
# summary
@test !isa(try summary(RawData()) catch ex ex end, Exception)
@test !isa(try summary(FFTData()) catch ex ex end, Exception)
@test !isa(try summary(CorrData()) catch ex ex end, Exception)
end
end
