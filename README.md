# SeisNoise.jl :sound: :earth_americas:
SeisNoise.jl is designed for fast and easy ambient noise cross-correlation on the CPU and GPU in Julia.

| **Documentation**                       | **Build Status**              | **Coverage** | **Chat**   |
|:---------------------------------------:|:-----------------------------------------:|:---------------------:|:---------------------:|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://tclements.github.io/SeisNoise.jl/latest) | [![Build Status](https://travis-ci.org/tclements/SeisNoise.jl.svg?branch=master)](https://travis-ci.org/tclements/SeisNoise.jl) |  [![Coverage Status](https://coveralls.io/repos/github/tclements/SeisNoise.jl/badge.svg?branch=master)](https://coveralls.io/github/tclements/SeisNoise.jl?branch=master) | [![](https://img.shields.io/badge/chat-on%20slack-yellow.svg)](https://slackinvite.julialang.org/) |


 ![Noise.jl Logo](/docs/src/assets/logo.png)

## Installation
You can install the latest version of SeisNoise using the Julia package manager (Press `]` to enter `pkg`).
From the Julia command prompt:

```julia
julia>]
(@v1.5) pkg> add SeisNoise
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("SeisNoise")
```

We recommend using the latest version of SeisNoise by updating with the Julia package manager:

```julia
(@v1.5) pkg> update SeisNoise
```

## Package Features

![flow](/docs/src/assets/SeisNoise-DataFlow.jpg)

  - Built upon [SeisIO](https://seisio.readthedocs.io/en/latest/) for easy and fast I/O.
  - Custom structures for storing Raw Data, Fourier Transforms of data, and cross-correlations
  - CPU/GPU compatible functions for cross-correlation.
  - Methods for [*dv/v* measurements](https://github.com/tclements/SeisDvv.jl).
  - Coming soon: Dispersion analysis.

## SeisNoise Cross-Correlation Example
Once you have installed the package you can type `using SeisNoise` to start
cross-correlating. SeisNoise uses a functional syntax to implement cross-correlation. For example

```Julia
using SeisNoise, SeisIO
fs = 40. # sampling frequency in Hz
freqmin,freqmax = 0.1,0.2 # min and max frequencies
cc_step, cc_len = 450, 1800 # corrleation step and length in S
maxlag = 60. # maximum lag time in correlation
s = "2019-02-03"
t = "2019-02-04"
S1 = get_data("FDSN","CI.SDD..BHZ",src="SCEDC",s=s,t=t)
S2 = get_data("FDSN","CI.PER..BHZ",src="SCEDC",s=s,t=t)
process_raw!(S1,fs)
process_raw!(S2,fs)
R = RawData.([S1,S2],cc_len,cc_step)
detrend!.(R)
taper!.(R)
bandpass!.(R,freqmin,freqmax,zerophase=true)
FFT = rfft.(R)
whiten!.(FFT,freqmin,freqmax)
C = correlate(FFT[1],FFT[2],maxlag)
clean_up!(C,freqmin,freqmax)
abs_max!(C)
corrplot(C)
```
will produce this figure:

![plot1](/docs/src/assets/xcorr-example.png)

## Cross-correlation on the GPU

SeisNoise can process data and compute cross-correlations on the GPU with CUDA. The [JuliaGPU](https://github.com/JuliaGPU) suite provides a high-level interface for CUDA programming through the CUDA.jl package. CUDA.jl provides an the `CuArray` type for storing data on the GPU. Data in SeisNoise structures (`R.x`, `F.fft`, and `C.corr` fields, for `RawData`, `FFTData`, and `CorrData`, respectively) can move between an `Array` on the CPU to a `CuArray` on the GPU using the `gpu` and `cpu` functions, as shown below.   

> :warning: Only **Nvidia** GPUs are suported at the moment. Hold in there for AMD/OpenCL support...

```julia
# create raw data and send to GPU
R = RawData(S1, cc_len, cc_step) |> gpu
R.x
72000×188 CUDA.CuArray{Float32,2,Nothing}

# send data back to the CPU
R = R |> cpu
R.x
72000×188 Array{Float32,2}
```

All basic processing remains the same on the GPU. Here is a complete cross-correlation routine on the GPU:

```julia
# send data to GPU
R1 = RawData(S1, cc_len, cc_step) |> gpu
R2 = RawData(S2, cc_len, cc_step) |> gpu
R = [R1,R2]

# preprocess on the GPU
detrend!.(R)
taper!.(R)
bandpass!.(R,freqmin,freqmax,zerophase=true)

# Real FFT on GPU
FFT = rfft.(R)
whiten!.(FFT,freqmin,freqmax)

# compute correlation and send to cpu
C = correlate(FFT[1],FFT[2],maxlag) |> cpu
```

### Routines Implemented on the GPU

![gpu times](/docs/src/assets/Fig2.jpg)

Processing times for a selection of routines on the GPU with Julia + GPU (white), Julia + CPU (black), and Python (grey). Currently these operations are implemented in SeisNoise on the GPU: 


- Preprocessing:
  - `detrend`,`demean`, `taper`, `onebit`, `smooth`
- Filtering:
  - `bandpass`, `bandstop`, `lowpass`, `highpass`
- Fourier Domain:
  - `whiten`, `rfft`, `irfft`
- Cross-correlation:
  - `correlate`, `cross-coherence`, `deconvolution`
- Post-processing:
  - `stack`, `filter`s, etc..

## Contributing
We welcome folks interested in contributing to SeisNoise. Please [open an issue](https://github.com/tclements/SeisNoise.jl/issues/new) to let us know about bug reports, new methods/code, and or feature requests/usage cases. If you would like to submit a pull request (PR), please include accompanying [tests](https://github.com/tclements/SeisNoise.jl/tree/master/test).
