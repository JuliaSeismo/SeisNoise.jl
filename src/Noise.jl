module Noise

using Dates, DataFrames, DSP, FFTW, JLD2, LinearAlgebra, Plots, SeisIO
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
include("compute_fft.jl")
include("correlation.jl")

end # module
