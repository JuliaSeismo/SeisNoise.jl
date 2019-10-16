push!(LOAD_PATH,"../src/")
using Documenter, SeisNoise

makedocs(
    modules = [SeisNoise],
    sitename = "SeisNoise.jl",
    # Uncomment below for local build
    # format = Documenter.HTML(prettyurls = false),
    authors = "Tim Clements",
    pages = Any[
        "Home" => "index.md",
        "Types" => "types.md",
        "Pre-Processing" => "preprocessing.md",
        "Filtering" => "filter.md",
        "Computing FFTs" => "fft.md",
        "Correlation" => "correlation.md",
        "Velocity Change" => "postprocessing.md",
        ],
)

deploydocs(
    repo = "github.com/tclements/SeisNoise.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
)
