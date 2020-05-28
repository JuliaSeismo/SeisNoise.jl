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
        "Workflow" => "types.md",
        "Using the GPU" => "using_gpus.md",
        "Pre-Processing" => "preprocessing.md",
        "Computing FFTs" => "fft.md",
        "Correlation" => "correlation.md",
        ],
)

deploydocs(
    repo = "github.com/tclements/SeisNoise.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
)
