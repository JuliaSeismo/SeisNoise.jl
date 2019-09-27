# SeisNoise.jl
SeisNoise.jl is designed for fast and easy ambient noise cross-correlation in Julia.

 [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://tclements.github.io/SeisNoise.jl/latest) [![Build Status](https://travis-ci.org/tclements/SeisNoise.jl.svg?branch=master)](https://travis-ci.org/tclements/SeisNoise.jl) [![Coverage Status](https://coveralls.io/repos/github/tclements/SeisNoise.jl/badge.svg?branch=master)](https://coveralls.io/github/tclements/SeisNoise.jl?branch=master)

 ![Noise.jl Logo](/docs/src/assets/logo.png)

## Installation
From the Julia command prompt:
1. Press `]` to enter `pkg`.
2. Type or copy: `add SeisNoise`
3. Press backspace to exit `pkg`.
4. Type or copy: `using SeisNoise`

## Package Features
  - Built upon [SeisIO](https://seisio.readthedocs.io/en/latest/) for easy and fast I/O.
  - Custom types for saving Fourier Transforms of data and cross-correlations
  - Array-based processing of raw data and cross-correlation.
  - Methods for *dv/v* measurements.
  - Coming soon: GPU support and dispersion analysis.

## Basic Cross-Correlation
Once you have installed the package you can type `using SeisNoise` to start
cross-correlating. For example

```Julia
using SeisNoise, SeisIO
fs = 10. # sampling frequency in Hz
freqmin,freqmax = 0.1,0.2 # minimum and maximum frequencies in Hz
cc_step, cc_len = 450, 1800 # corrleation step and length in S
maxlag = 60. # maximum lag time in correlation
smoothing_half_win = 3
S1 = get_data("FDSN","CI.SDD..BHZ",src="SCEDC",s="2019-02-03",t="2019-02-04")
S2 = get_data("FDSN","CI.PER..BHZ",src="SCEDC",s="2019-02-03",t="2019-02-04")
FFT1 = compute_fft(S1,freqmin, freqmax, fs, cc_step, cc_len)
FFT2 = compute_fft(S2,freqmin, freqmax, fs, cc_step, cc_len)
coherence!(FFT1,smoothing_half_win)
coherence!(FFT2,smoothing_half_win)
C = compute_cc(FFT1,FFT2,maxlag,corr_type="coherence")
clean_up!(C,freqmin,freqmax)
abs_max!(C)
corrplot(C)
```
will produce this figure:

![plot1](/docs/src/assets/xcorr-example.png)
