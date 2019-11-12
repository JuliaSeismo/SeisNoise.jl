__precompile__()
module SeisNoise

using Dates, DataFrames, DSP, FFTW, Glob, JLD2, LinearAlgebra, SeisIO
using Statistics, StatsBase, Interpolations, GLM, Plots, Distributed, LightXML

# import types first
include("Types/RawData.jl")
include("Types/FFTData.jl")
include("Types/CorrData.jl")
include("Types/show.jl")
include("Types/InputParams.jl")

# import pre and post processing tools
include("distance.jl")
include("io.jl")
include("ArrayFuncs.jl")
include("tools.jl")
include("slicing.jl")
include("filter.jl")
include("downsample.jl")
include("availability.jl")
include("phase_shift.jl")
include("stacking.jl")

# import  routines for doin' stuff
include("compute_fft.jl")
include("correlation.jl")
include("rotation.jl")
include("VelocityChange/MWCS.jl")
include("VelocityChange/Stretching.jl")
include("VelocityChange/Wavelets.jl")
include("Plotting/plotting.jl")

end # module
