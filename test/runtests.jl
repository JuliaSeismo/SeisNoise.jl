using SeisNoise, SeisIO, Dates, FFTW, Statistics
using Test

include("test_NoiseData.jl")
include("test_show.jl")
include("test_ArrayFuncs.jl")
include("test_distance.jl")
# include("test_io.jl")
include("test_tools.jl")
include("test_slicing.jl")
include("test_filter.jl")
include("test_phase_shift.jl")
# include("test_stacking.jl")
#
# # import  routines for doin' stuff
# include("test_compute_fft.jl")
# include("test_correlation.jl")
# include("test_rotation.jl")
# include("test_plotting.jl")
#
# # out with the old
# include("test_deprecated.jl")
