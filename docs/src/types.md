# SeisNoise Types

`SeisNoise` is designed around array-based cross-correlation. `SeisNoise` uses three custom data structures: `RawData` for storing windows ambient noise, `FFTData` for
storing Fourier transforms of ambient noise, and `CorrData` for storing
cross-correlations. SeisNoise types extend `SeisIO`'s `SeisChannel` type
for 2D-array ambient noise processing. SeisNoise's modular functions work on `RawData`, `FFTData`, and `CorrData` objects *in-place* through Julia's
[multiple-dispatch](https://en.wikipedia.org/wiki/Multiple_dispatch). Functions in
SeisNoise that end in ! (e.g. `onebit!(R)`) by convention modify their arguments,
while functions that do not end in ! (e.g. `onebit(R)`) allocate new arrays or
objects.

!!! note

      By convention in Julia, array data is stored *column-wise*. Noise windows in
      `RawData` objects, ffts in `FFTData` objects, and cross-corrleations in
      `CorrData` objects are all stored column-wise. To access data one or more time windows, use column-based indexing, e.g. `R[:,1]`, `F[:,2]`, or `C[:,3:5]`.

## `RawData` - Objects for windowed raw data
`RawData` objects store windows of ambient noise. `RawData` have very similar structure to `SeisChannel` objects, except with added fields for:
- `cc_len`: the cross-correlation window length (in seconds).
- `cc_step`: the step length between successive correlation windows (in seconds).
- `freqmin`: minimum frequency `RawData` has been filtered (in Hz).
- `freqmax`: maximum frequency `RawData` has been filtered (in Hz).
- `time_norm`: One-bit whitening or phase-whitening applied.
- `t`: Starttime of each noise window (number of seconds since the unix epoch 1970-01-01T00:00:00 as
  a Float64).
- `x`: Raw data stored in columns (2d-array).


To create an empty `RawData` object, use the `RawData()` function. Here is an empty `RawData` object in the REPL:

```julia
julia> using SeisNoise
julia> RawData()
RawData with 0 windows
      NAME: …
        ID: …
       LOC: 0.0 N, 0.0 E, 0.0 m
        FS: 0.0
      GAIN: 1.0
   FREQMIN: 0.0
   FREQMAX: 0.0
    CC_LEN: …
   CC_STEP: …
 TIME_NORM: …
      RESP: a0 1.0, f0 1.0, 0z, 0p
      MISC: …
     NOTES: …
         T: …
         X: …
```

A more natural way to create a `RawData`object is to convert an existing `SeisChannel`
or `SeisData` object into a `RawData` object. This allows one to do preprocessing
that is more natural for single component data (e.g. instrument response removal) before
doing array-based processing. Here is an example of converting a `SeisData` object
to a `RawData` object.

!!! note

      `cc_len` and `cc_step` *must* be included to create a `RawData` object.

```julia
julia> using SeisIO, SeisNoise

julia> S = get_data("FDSN", "TA.M22K..BHZ", src="IRIS",s="2019-01-01",t=600)
SeisData with 1 channels (1 shown)
    ID: TA.M22K..BHZ                       
  NAME: Willow, AK, USA                    
   LOC: 61.7531 N, -150.12 E, 57.0 m       
    FS: 40.0                               
  GAIN: 5.03883e8                          
  RESP: a0 8.32666e17, f0 0.2, 6z, 11p     
 UNITS: m/s                                
   SRC: http://service.iris.edu/fdsnws/da…
  MISC: 4 entries                          
 NOTES: 0 entries                          
     T: 2019-01-01T00:00:00.000 (0 gaps)   
     X: +8.540e+02                         
        +6.230e+02                         
            ...                            
        +1.081e+03                         
        (nx = 24001)                       
     C: 0 open, 0 total

julia> cc_len, cc_step = 100,50 # specify window length and step (in seconds)
(100, 50)

julia> R = RawData(S,cc_len,cc_step)
RawData with 11 windows
      NAME: "Willow, AK, USA"                  
        ID: "TA.M22K..BHZ"                     
       LOC: 61.7531 N, -150.12 E, 57.0 m
        FS: 40.0
      GAIN: 5.03883e8
   FREQMIN: 0.01
   FREQMAX: 20.0
    CC_LEN: 100                                
   CC_STEP: 50                                 
 TIME_NORM: false                              
      RESP: a0 8.32666e17, f0 0.2, 6z, 11p
      MISC: 4 entries                          
     NOTES: 2 entries                          
         T: 2019-01-01T00:00:00.000            …
         X: 4000×11 Array{Float32,2}     

```

This created a `RawData` object, `R`, 11 windows, each 100 seconds long and sampled at `40` Hz (4000 points). The `freqmin` is automatically set to `0.01` Hz because this is the `1 / cc_len`, whereas `freqmax` is set to the Nyquist frequency `R.fs/2`. To access the noise data stored in `R` just do `R.x`

```julia
julia> R.x
4000×11 Array{Float32,2}:
  854.0  1645.0   790.0  976.0   959.0  …   921.0  1148.0   718.0
  623.0   292.0   762.0  790.0  1065.0      828.0  1242.0   596.0
  530.0  1398.0   408.0  649.0  1033.0      -32.0  1263.0   678.0
  684.0   970.0   546.0  793.0   984.0     1399.0  1324.0   676.0
  810.0   273.0  1066.0  943.0   992.0      561.0  1340.0   655.0
  834.0   554.0  1385.0  866.0   955.0  …    62.0  1291.0   729.0
    ⋮                                   ⋱                     ⋮  
  647.0   810.0   801.0  380.0  1683.0  …   830.0   815.0  1223.0
  639.0   700.0   880.0  285.0  1284.0      709.0   798.0  1149.0
  860.0   789.0   955.0  133.0  1427.0      844.0   913.0   998.0
 1052.0   875.0   927.0  115.0  1623.0      857.0   944.0   970.0
  872.0   925.0   883.0    2.0  1141.0      810.0   883.0  1061.0
```

See the Pre-Processing and Filtering pages for
functions that work on `RawData` objects.  


## `FFTData` - Objects for Fourier transforms (FFTs)

Cross-correlation in `SeisNoise` is done in the frequency domain:

```math
C_{AB}(ω) = u^*_A(ω) u_B(ω)
```

where ``C_{AB}(ω)`` is the element-wise multiplication between
``u^*_A(ω)``, the complex conjugate of the Fourier transform of ambient noise time series ``A``, and
``u_B(ω)``, the Fourier transform of ambient noise
time series ``B``. This exploits the **O(nlogn)** complexity of the Fast-Fourier Transforms. Since ambient noise data
is real valued, `SeisNoise` uses real Fast-Fourier Transforms (`rfft`), which offer a 2-3x
speed over a general `fft`.

The `FFTData` structure in `SeisNoise` stores
Fourier transforms (``u(ω)``) of ambient noise. `FFTData` allow users to apply smoothing
operations, such as whitening, coherence, or deconvolution, in-place before
cross-corrleating. To create an empty `FFTData` object, use the `FFTData()`
function.

```julia
julia> using SeisNoise
julia> FFTData()
FFTData with 0 ffts
      NAME: …
        ID: …
       LOC: 0.0 N, 0.0 E, 0.0 m
        FS: 0.0
      GAIN: 1.0
   FREQMIN: 0.0
   FREQMAX: 0.0
    CC_LEN: …
   CC_STEP: …
  WHITENED: …
 TIME_NORM: …
      RESP: c = 0.0, 0 zeros, 0 poles
      MISC: …
     NOTES: …
         T: …
       FFT: …

```
The only difference between an `FFTData` and a `RawData` object is the addition
of the `whitened` parameter and the swap of the `x` data field to the `fft` data
field.

`FFTData` more naturally flow from `RawData` input to the
 `compute_fft` function.

 ```julia
 julia> F = compute_fft(R)
FFTData with 11 ffts
      NAME: "TA.M22K..BHZ"                     
        ID: "2019-01-01"                       
       LOC: 61.7531 N, -150.12 E, 57.0 m
        FS: 40.0
      GAIN: 5.03883e8
   FREQMIN: 0.01
   FREQMAX: 20.0
    CC_LEN: 100                                
   CC_STEP: 50                                 
  WHITENED: false                              
 TIME_NORM: false                              
      RESP: a0 8.32666e17, f0 0.2, 6z, 11p
      MISC: 4 entries                          
     NOTES: 2 entries                          
         T: 2019-01-01T00:00:00.000            …
       FFT: 2001×11 Array{Complex{Float32},2}  
```

To access the Fourier transform stored in `F` just do `F.fft`

```julia
julia> F.fft
2001×11 Array{Complex{Float32},2}:
 3.18286e6+0.0im      3.18685e6+0.0im       …  3.19561e6+0.0im    
  -3205.34-7445.99im    12079.2+1573.37im        11842.9-4106.75im
  -10263.9+5427.09im    -1056.6+3629.68im        5421.51+2751.05im
  -6255.81-25112.9im    21741.7+9154.63im        37171.6+10855.5im
  -6164.36-7888.87im    27625.1-2106.65im        10308.3-24606.1im
   2470.13+8848.62im    13407.9+5940.08im   …   -10548.2+1305.7im
          ⋮                                 ⋱           ⋮         
   133.262-18.7686im    1204.43+18.5029im       -135.982+7.42407im
   122.472-7.24219im    1179.84+1.39209im       -133.742-5.74707im
   131.036+7.63867im    1199.29+0.913574im      -117.602+4.26135im
   123.142+4.82568im    1202.13+1.8291im        -112.404+11.0547im
     126.0+0.0im         1179.0+0.0im       …     -148.0+0.0im    
```


## `CorrData` - Objects for ambient noise cross-correlations

As stated in the previous section, cross-correlation in the frequency domain is
an element-wise product of Fourier-transforms of ambient noise:

```math
C_{AB}(ω) = u^*_A(ω) u_B(ω)
```

To transform a frequency domain cross-correlation into the time domain, we take
the inverse real Fourier transform of ``C_{AB}(ω)``:

```math
C(τ)_{AB} = \mathfrak{F}^{-1} \left(C_{AB}(ω)\right)
```

where ``τ`` is the lag time. Computing a cross-correlation thus necessitates two
Fourier transforms, one element-wise multiplication (of complex numbers), and an
inverse Fourier transform.


The `CorrData` structure in `SeisNoise` stores time-domain cross-correlations.
To create an empty `CorrData` object, use the `CorrData()` function.

```julia
julia> CorrData()
CorrData with 0 Corrs
      NAME: …
        ID: …
       LOC: 0.0 N, 0.0 E, 0.0 m
      COMP: …
   ROTATED: …
 CORR_TYPE: …
        FS: 0.0
      GAIN: 1.0
   FREQMIN: 0.0
   FREQMAX: 0.0
    CC_LEN: …
   CC_STEP: …
  WHITENED: …
 TIME_NORM: …
      RESP: a0 1.0, f0 1.0, 0z, 0p
      MISC: …
     NOTES: …
      DIST: 0.0
       AZI: 0.0
       BAZ: 0.0
    MAXLAG: 0.0
         T: …
      CORR: …
```
`CorrData` have a few additional parameters:
- `comp`: The component of the cross-correlations, e.g. EN, RT, ZZ, etc...
- `rotated`: True/false if correlation has been rotated.
- `corr_type`: The type of correlation, e.g. "cross-correlation", "coherence", or "deconvolution".
- `dist`: Distance between station 1 and station 2 (in Km).
- `azi`: Azimuth from station 1 and station 2 (in degrees).
- `baz`: Back azimuth between station 1 and station 2 (in degrees).
- `maxlag`: Maximum lag time (in seconds) in cross-correlations.
- `corr`: Time-domain cross-correlations stored in columns.

`CorrData` more naturally flow from `FFTData` input to the `compute_cc` function. In this example,
we are computing an "autocorrleation" by inputting the same an `FFTData` twice:

```julia
julia> maxlag = 20.
20.0

julia> C = compute_cc(F,F,maxlag,corr_type="cross-correlation")
CorrData with 11 Corrs
      NAME: "TA.M22K..BHZ.TA.M22K..BHZ"        
        ID: "2019-01-01"                       
       LOC: 61.7531 N, -150.12 E, 57.0 m
      COMP: "ZZ"                               
   ROTATED: false                              
 CORR_TYPE: "cross-correlation"                
        FS: 40.0
      GAIN: 5.03883e8
   FREQMIN: 0.01
   FREQMAX: 20.0
    CC_LEN: 100                                
   CC_STEP: 50                                 
  WHITENED: false                              
 TIME_NORM: false                              
      RESP: a0 8.32666e17, f0 0.2, 6z, 11p
      MISC: 4 entries                          
     NOTES: 2 entries                          
      DIST: 0.0
       AZI: 0.0
       BAZ: 0.0
    MAXLAG: 20.0
         T: 2019-01-01T00:00:00.000            …
      CORR: 1601×11 Array{Float32,2}   
```

To access the autocorrelations stored in `C` just do `C.corr`:

```julia
julia> C.corr
1601×11 Array{Float32,2}:
 2.56154e9  2.5478e9   2.51531e9  …  2.53951e9  2.55251e9  2.51616e9
 2.54715e9  2.52733e9  2.51668e9     2.54039e9  2.55533e9  2.51614e9
 2.53042e9  2.54303e9  2.52005e9     2.53713e9  2.55813e9  2.51694e9
 2.57688e9  2.57981e9  2.52313e9     2.53812e9  2.55982e9  2.52244e9
 2.617e9    2.57287e9  2.52175e9     2.53825e9  2.56081e9  2.52769e9
 2.51533e9  2.50786e9  2.5218e9   …  2.53502e9  2.56186e9  2.52965e9
 2.5034e9   2.52427e9  2.52484e9     2.53588e9  2.56324e9  2.53052e9
 2.5998e9   2.578e9    2.52705e9     2.53676e9  2.56418e9  2.53246e9
 ⋮                                ⋱                        ⋮        
 2.56922e9  2.56349e9  2.52909e9     2.53253e9  2.56514e9  2.53592e9
 2.5998e9   2.578e9    2.52705e9  …  2.53676e9  2.56418e9  2.53246e9
 2.5034e9   2.52427e9  2.52484e9     2.53588e9  2.56324e9  2.53052e9
 2.51533e9  2.50786e9  2.5218e9      2.53502e9  2.56187e9  2.52965e9
 2.617e9    2.57287e9  2.52175e9     2.53825e9  2.56081e9  2.52769e9
 2.57688e9  2.57981e9  2.52313e9     2.53812e9  2.55982e9  2.52244e9
 2.53042e9  2.54303e9  2.52005e9  …  2.53713e9  2.55813e9  2.51694e9
```

### Type Documentation
```@docs
RawData
FFTData
CorrData
```
