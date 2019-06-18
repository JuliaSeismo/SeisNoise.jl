push!(LOAD_PATH,"../src/")
using Documenter, Noise

makedocs(
    modules = [Noise],
    sitename = "Noise.jl",
    authors = "Tim Clements",
    pages = Any[
        "Home" => "index.md",
        "Pre-Processing" => "preprocessing.md",
        "Array Functions" => "arrayfuncs.md",
        "Filters" => "filter.md",
        "Computing FFTs" => "fft.md",
        "Correlation" => "correlation.md",
        "Velocity Change" => "postprocessing.md",
        "Types" => Any[
         "InputParams" => "Types/inputparams.md",
         "FFTData" => "Types/fftdata.md",
         "CorrData" => "Types/corrdata.md",
        ],
        ],
)

deploydocs(
    repo = "github.com/tclements/Noise.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
)
