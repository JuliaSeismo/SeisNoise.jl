# `Computing Correlations` - methods for computing correlations from FFTs.

Cross-correlation in the frequency domain is element-wise multiplication between
``u_A(ω)``, the Fourier transform of ambient noise time series ``A``, and
``u^*_B(ω)``, the complex conjugate of the Fourier transform of ambient noise
time series ``B``. Options for cross-correlation in SeisNoise.jl include

cross-correlation:

```math
C_{AB}(ω) = u_A(ω) u^∗_B(ω)
```

cross-coherency:

```math
C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_A(\omega) ∣ \, ∣ u_B(\omega) ∣}
```

and deconvolution:

```math
C_{AB}(\omega) = \frac{u_A(\omega) u^∗_B(\omega)}{\mid u_B(\omega) {\mid}^2}
```

The cross-correlation in the time domain is just the inverse real Fourier transform of ``C_{AB}(ω)``:

```math
C(τ)_{AB} = \mathfrak{F}^{-1} \left(C_{AB}(ω)\right)
```

where ``τ`` is the lag time.

## Computing Correlations

The `compute_cc` function provides the typical workflow for computing correlations
in SeisNoise.jl. The necessary inputs to `compute_cc` are the maximum lag time, `max_lag`,
in seconds to save, e.g. 200 seconds, the type of correlation( e.g. `cross-correlate`,
`coherence`, or `deconv`), and the number of points to smooth the spectrum of
``u_A(ω)`` or ``u_B(ω)`` if using the coherence or deconvolution.

```julia
using SeisNoise, SeisIO
julia> fs = 40. # sampling frequency in Hz
julia> freqmin,freqmax = 0.1,0.2 # minimum and maximum frequencies in Hz
julia> cc_step, cc_len = 450, 1800 # corrleation step and length in S
julia> maxlag = 80. # maximum lag time in correlation
julia> S1 = get_data("IRIS","TA.V04C..BHZ",s="2006-02-01",t="2006-02-02")
julia> S2 = get_data("IRIS","TA.V05C..BHZ",s="2006-02-01",t="2006-02-02")
julia> FFT1 = compute_fft(S1,freqmin, freqmax, fs, cc_step, cc_len,
                  time_norm=false,to_whiten=false)
julia> FFT2 = compute_fft(S2,freqmin, freqmax, fs, cc_step, cc_len,
                  time_norm=false,to_whiten=false)
julia> C = compute_cc(FFT1,FFT2,maxlag,corr_type="coherence")
CorrData with 188 Corrs
      NAME: "TA.V04C..BHZ.TA.V05C..BHZ"        
        ID: "2006-02-01"                       
       LOC: 0.0 N, 0.0 E, 0.0 m
      COMP: "ZZ"                               
   ROTATED: false                              
 CORR_TYPE: "coherence"                        
        FS: 40.0
      GAIN: 1.0
   FREQMIN: 0.1
   FREQMAX: 0.2
    CC_LEN: 1800                               
   CC_STEP: 450                                
  WHITENED: false                              
 TIME_NORM: false                              
      RESP: c = 0.0, 0 zeros, 0 poles
      MISC: 0 entries                          
     NOTES: 2 entries                          
    MAXLAG: 80.0
         T: 2006-02-01T00:07:30.000            …
      CORR: 6401×188 Array{Float32,2}
```

## Saving/Loading Correlations

`CorrData` objects can be saved to disk in the native Julia [JLD2](https://github.com/JuliaIO/JLD2.jl)
format using the `save_corr` function.

```julia
julia> OUTDIR = "~/TEST/CORR/"
julia> save_corr(C,OUTDIR)
```

`CorrData` are stored in groups by component (e.g. ZZ or RZ), then by date
(in yyyy-mm-dd format) in JLD2. By default, JLD2 files are saved to
/PATH/TO/OUTDIR/NET1.STA1.CHAN1.NET2.STA2.CHAN2.jld2.

```julia
file = jldopen("~/TEST/CORR/TA.V04C.BHZ.TA.V05C.BHZ.jld2","r")
JLDFile ~/TEST/CORR/TA.V04C.BHZ.TA.V05C.BHZ.jld2 (read-only)
 └─📂 ZZ
    └─🔢 2006-02-01
```

To read an `CorrData` on disk, use the `load_corr` function:

```julia
julia> C = load_corr("~/TEST/CORR/TA.V04C.BHZ.TA.V05C.BHZ.jld2","ZZ")
CorrData with 188 Corrs
      NAME: "TA.V04C..BHZ.TA.V05C..BHZ"        
        ID: "2006-02-01"                       
       LOC: 0.0 N, 0.0 E, 0.0 m
      COMP: "ZZ"                               
   ROTATED: false                              
 CORR_TYPE: "coherence"                        
        FS: 40.0
      GAIN: 1.0
   FREQMIN: 0.1
   FREQMAX: 0.2
    CC_LEN: 1800                               
   CC_STEP: 450                                
  WHITENED: false                              
 TIME_NORM: false                              
      RESP: c = 0.0, 0 zeros, 0 poles
      MISC: 0 entries                          
     NOTES: 2 entries                          
    MAXLAG: 80.0
         T: 2006-02-01T00:07:30.000            …
      CORR: 6401×188 Array{Float32,2}  
```

```@docs
clean_up!
correlate
compute_cc
save_corr
load_corr
stack!
```
