__precompile__()
module Noise

using Dates, DataFrames, DSP, FFTW, Glob, JLD2, LinearAlgebra, SeisIO
# include("ArrayFuncs.jl")
include("tools.jl")
include("slicing.jl")
include("filter.jl")
include("downsample.jl")
include("availability.jl")
include("phase_shift.jl")
include("Types/FFTData.jl")
include("Types/CorrData.jl")
include("Types/show.jl")
include("Types/InputParams.jl")
include("compute_fft.jl")
include("correlation.jl")

end # module
