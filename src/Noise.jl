module Noise

using Dates, DSP, FFTW, LinearAlgebra, Plots, PyCall, Statistics
include("tools.jl")
include("filter.jl")
include("correlate.jl")
include("availability.jl")
include("compute_cc.jl")

end # module
