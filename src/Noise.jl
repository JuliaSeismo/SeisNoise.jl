__precompile__()
module Noise

using Dates, DataFrames, DSP, FFTW, Glob, JLD2, LinearAlgebra, Printf, SeisIO
using  Statistics, Interpolations, GLM, Plots

# import types first
include("Types/FFTData.jl")
include("Types/CorrData.jl")
include("Types/show.jl")
include("Types/InputParams.jl")

# import pre and post processing tools
include("ArrayFuncs.jl")
include("tools.jl")
include("slicing.jl")
include("filter.jl")
include("downsample.jl")
include("availability.jl")
include("phase_shift.jl")

# import  routines for doin' stuff
include("compute_fft.jl")
include("correlation.jl")
include("VelocityChange/MWCS.jl")
include("VelocityChange/Stretching.jl")
include("Plotting/plotting.jl")

end # module
