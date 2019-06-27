# Noise.jl

*Ambient Noise Cross-Correlation in Julia.*

Noise.jl provides routines for quickly and efficiently implementing [seismic interferometry](https://en.wikipedia.org/wiki/Seismic_interferometry).

!!! note

    Much of the I/O and preprocessing in Noise.jl uses the [SeisIO.jl](https://github.com/jpjones76/SeisIO.jl) package. Please read through the [SeisIO Documentation](https://seisio.readthedocs.io/en/latest/) to get familiar with seismic data processing in Julia.

## Installation
From the Julia command prompt:
1. Press `]` to enter `pkg`.
2. Type or copy: `add https://github.com/tclements/Noise.jl; build; precompile`
3. Press backspace to exit `pkg`.
4. Type or copy: `using Noise`

## Using Noise.jl

Noise.jl was designed to be as easy to use in the REPL as on an HPC cluster. If
you want to get started processing data, head over to the tutorial or parallel example.
This documentation provides a reference to all the underlying function for cross-correlation.
We encourage you to develop your own workflow using Noise's core functionality.

![plot1](assets/CI-moveout.png)
