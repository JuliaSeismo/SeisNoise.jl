module SeisNoise

using Dates, DSP, FFTW, Glob, JLD2, LinearAlgebra, SeisBase, SeisBase.Nodal
using Statistics, StatsBase, CUDA, Adapt, GPUArrays, RecipesBase

# check use of cuda
const use_cuda = Ref(false)
if !CUDA.functional()
  else
    use_cuda[] = true
end

# import types first
include("Types/NoiseData.jl")
include("Types/NodalData.jl")
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
include("nodalcorrelation.jl")
include("rotation.jl")
include("Plotting/plotting.jl")

# out with the old
include("deprecated.jl")

end # module
