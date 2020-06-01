# Computing FFTs

All correlation in SeisNoise.jl is done in the frequency domain, which can be represented
by:

```math
C_{AB}(ω) = u_A(ω) u^∗_B(ω)
```

Where ``u_A(ω)`` is the Fourier transform of ambient noise time series ``A``,
``u^*_B(ω)`` is the complex conjugate of the Fourier transform of ambient noise
time series ``B``, and ``C_{AB}(ω)`` is the cross-correlation of ``A`` and ``B``
in the frequency domain. For time and memory efficiency, the real Fourier transform
(`rfft`) is used as opposed to the regular Fourier transform (`fft`). This gives
a speedup of about a factor of 3. A typical workflow for computing `fft`'s can include:

- Real Fourier Transform
- Spectral whitening (removing spectral amplitude information)
- Hilbert Transform (for phase-correlation)

## Creating FFTData Structures

The `rfft` function provides the typical workflow for computing `fft`s in SeisNoise.

```julia
using SeisNoise, SeisIO
S = get_data("IRIS","TA.V04C..BHZ",s="2006-02-01T00:00:00",t="2006-02-01T01:00:00")
cc_step, cc_len = 100., 100.
R = RawData(S,cc_len,cc_step)
F = rfft(R)
FFTData with 35 ffts
      NAME: "TA.V04C..BHZ"                     
        ID: "2006-02-01"                       
       LOC: 0.0 N, 0.0 E, 0.0 m
        FS: 40.0
      GAIN: 1.0
   FREQMIN: 0.01
   FREQMAX: 20.0
    CC_LEN: 100.0
   CC_STEP: 100.0
  WHITENED: false                              
 TIME_NORM: ""                                 
      RESP: a0 1.0, f0 1.0, 0z, 0p
      MISC: 0 entries                          
     NOTES: 2 entries                          
         T: 2006-02-01T00:01:40                …
       FFT: 2001×35 Array{Complex{Float32},2}  
```

Underneath the hood, `rfft` applies a real Fourier transform to the `.x` field of a
`RawData` structure, then allocates a new `FFTData` structure with the Fourier transform
data stored in the `.fft` field.

## Whitening FFTData

SeisNoise provides three methods for normalizing the spectrum of an FFTData structure.
The `whiten!` method, sets the complex amplitude of frequencies between `freqmin`
and `freqmax` to 1. This preserves only the phase component of the signal.  

```julia
freqmin, freqmax = 10., 20.
whiten!(F,freqmin,freqmax)
```

The `coherence` method smooths an the spectrum of an `FFTData` by the smoothed
absolute value of itself.

```julia
coherence!(F,20)
```

The `deconvolution` method smooths an the spectrum of an `FFTData` by the smoothed
absolute value squared of itself.

```julia
deconvolution!(F,20)
```

## Saving/Loading FFTs

`FFTData` objects can be saved to disk in the native Julia [JLD2](https://github.com/JuliaIO/JLD2.jl)
format using the `save_fft` function.

```julia
julia> OUTDIR = "~/TEST/FFT/"
julia> save_fft(F,OUTDIR)
```

`FFTData` are stored in groups by channel (e.g. BHZ or HHZ), then by date
(in yyyy-mm-dd format) in JLD2. By default, JLD2 files are saved to
/PATH/TO/OUTDIR/NET.STA.CHAN.jld2.

```julia
file = jldopen("~/TEST/FFT/TA.V04C.BHZ.jld2","r")
JLDFile ~TEST/FFT/TA.V04C.BHZ.jld2 (read-only)
 └─📂 BHZ
    └─🔢 2006-02-01
```

To read an `FFTData` on disk, use the `load_fft` function:

```julia
julia> F = load_fft("~/TEST/FFT/TA.V04C.BHZ.jld2","BHZ")
FFTData with 35 ffts
      NAME: "TA.V04C..BHZ"                     
        ID: "2006-02-01"                       
       LOC: 0.0 N, 0.0 E, 0.0 m
        FS: 20.0
      GAIN: 1.0
   FREQMIN: 0.05
   FREQMAX: 5.0
    CC_LEN: 100                                
   CC_STEP: 100                                
  WHITENED: false                              
 TIME_NORM: false                              
      RESP: c = 0.0, 0 zeros, 0 poles
      MISC: 0 entries                          
     NOTES: 2 entries                          
         T: 2006-02-01T00:00:00.000            …
       FFT: 2001×35 Array{Complex{Float32},2}
```

Note that it is necessary to specify the channel when using `load_fft`.

```@docs
whiten!
SeisNoise.rfft
coherence!
deconvolution!
save_fft
load_fft
```
