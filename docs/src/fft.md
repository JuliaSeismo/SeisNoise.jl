# `Computing FFTs` - methods for computing FFTs raw noise data.

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

- Spectral whitening (removing spectral amplitude information)
- One-bit normalization
- Phase normalization
- Real Fourier Transform

## Computing FFTs

The `compute_fft` function provides the typical workflow for computing `fft`s.

```julia
julia> using SeisNoise, SeisIO
julia> S = get_data("IRIS","TA.V04C..BHZ",s="2006-02-01T00:00:00",t="2006-02-01T01:00:00")
julia> freqmin, freqmax = 1.,10.
julia> fs = 20.
julia> cc_step, cc_len = 100, 100
julia> F = compute_fft(S,freqmin,freqmax,fs,cc_step,cc_len,time_norm=false,
                       to_whiten=false)

FFTData with 36 ffts
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
       FFT: 1001×36 Array{Complex{Float32},2}  
```

Underneath the hood, `compute_fft` applies pre-processing with `merge`, `ungap`,
and `process_raw`. The `SeisData` object is then transformed into an Array `A` of
sliding windows. The `process_fft` then applies spectral or time-domain normalization
and returns `FFT`, the Fourier transform of `A`. An `FFTData` object is then created
from the parameters of `S` and `FFT`.

```julia                        
merge!(S)
ungap!(S)
process_raw!(S,fs)
A, starts, ends = slide(S[1], cc_len, cc_step)
FFT = process_fft(A, freqmin, freqmax, fs, time_norm=time_norm,to_whiten=to_whiten)
F = FFTData(S[1].id, Dates.format(u2d(starts[1]),"Y-mm-dd"),
                       S[1].loc, S[1].fs, S[1].gain, freqmin, freqmax,
                       cc_len, cc_step, to_whiten, time_norm, S[1].resp,
                       S[1].misc, S[1].notes, starts, FFT)
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
FFTData with 36 ffts
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
       FFT: 1001×36 Array{Complex{Float32},2}
```

Note that it is necessary to specify the channel when using `load_fft`.

```@docs
whiten
process_fft
compute_fft
save_fft
load_fft
```
