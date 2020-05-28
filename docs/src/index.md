# SeisNoise.jl

*Ambient Noise Cross-Correlation in Julia.*

SeisNoise.jl provides routines for quickly and efficiently implementing [seismic interferometry](https://en.wikipedia.org/wiki/Seismic_interferometry).

!!! note

    Much of the I/O and preprocessing in SeisNoise.jl uses the [SeisIO.jl](https://github.com/jpjones76/SeisIO.jl) package. Please read through the [SeisIO Documentation](https://seisio.readthedocs.io/en/latest/) to get familiar with seismic data processing in Julia.



## Installation
You can install the latest version of SeisNoise using the Julia package manager (Press `]` to enter `pkg`).
From the Julia command prompt:

```julia
julia>]
(@v1.4) pkg> add SeisNoise
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("SeisNoise")
```

We recommend using the latest version of SeisNoise by updating with the Julia package manager:

```julia
(@v1.4) pkg> update SeisNoise
```

## Package Features
- Built upon [SeisIO](https://seisio.readthedocs.io/en/latest/) for easy and fast I/O.
- Custom structures for storing Raw Data, Fourier Transforms of data, and cross-correlations
- CPU/GPU compatible functions for cross-correlation.
- Methods for [*dv/v* measurements](https://github.com/tclements/SeisDvv.jl).
- Coming soon: Dispersion analysis.

## Getting Started
SeisNoise.jl was designed to be as easy to use in the REPL as on an HPC cluster. This documentation provides a reference to all the underlying function for cross-correlation.
We encourage you to develop your own workflow using SeisNoise's core functionality.

![plot1](assets/CI-moveout.png)
