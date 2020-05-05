module SeisNoise

using Dates, DataFrames, DSP, FFTW, Glob, JLD2, LinearAlgebra, SeisIO
using Statistics, StatsBase, Plots, Distributed
using Distributed, CuArrays, Adapt, CUDAnative, GPUArrays

# check use of cuda
const use_cuda = Ref(false)
if !CuArrays.functional()
  else
    use_cuda[] = true
end

# import types first
include("Types/NoiseData.jl")
include("Types/show.jl")

# import pre and post processing tools
include("distance.jl")
include("io.jl")
include("ArrayFuncs.jl")
include("tools.jl")
include("slicing.jl")
include("filter.jl")
include("phase_shift.jl")
include("stacking.jl")

# import  routines for doin' stuff
include("compute_fft.jl")
include("correlation.jl")
include("rotation.jl")
include("Plotting/plotting.jl")

# out with the old
include("deprecated.jl")

end # module
