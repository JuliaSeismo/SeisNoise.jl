module Noise

using Dates, DataFrames, DSP, FFTW, JLD2, LinearAlgebra, Plots, SeisIO
include("ArrayFuncs.jl")
include("tools.jl")
include("slicing.jl")
include("filter.jl")
include("downsample.jl")
include("correlate.jl")
include("availability.jl")
include("phase_shift.jl")
include("FFTData.jl")
include("compute_fft.jl")

end # module
