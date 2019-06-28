`FFTData` - Objects for holding Fourier transforms (FFTs)

`SeisNoise` is designed around array-based cross-correlation. `SeisNoise` uses a custom
structure `FFTData` for holding Fourier transforms of ambient noise. To create
an empty `FFTData` object, use the `FFTData()` function. `FFTData` are created
with the `compute_fft` function.

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

A non-empty `FFTData` looks like this:
```julia
julia> julia> FFT
FFTData with 188 ffts
      NAME: "TA.V04C..BHZ"                     
        ID: "2006-02-01"                       
       LOC: 0.0 N, 0.0 E, 0.0 m
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
         T: 2006-02-01T00:07:30.000            …
       FFT: 36001×188 Array{Complex{Float32}…     
```

To access the Fourier transformed stored in `FFT` just do `FFT.fft`

```julia
julia> FFT.fft
36001×188 Array{Complex{Float32},2}:
 1.55183e-11+0.0im          1.59162e-12+0.0im          …  -9.66338e-12+0.0im        
 -7.41018e-7-7.1289e-9im    -4.28328e-7-6.37775e-9im        4.01078e-7-2.87796e-9im
  -7.4106e-7-1.42547e-8im   -4.28291e-7-1.27635e-8im         4.0114e-7-5.75133e-9im
 -7.41128e-7-2.13963e-8im   -4.28259e-7-1.91511e-8im        4.01203e-7-8.65502e-9im
 -7.41227e-7-2.85335e-8im   -4.28155e-7-2.55307e-8im        4.01347e-7-1.15261e-8im
 -7.41355e-7-3.5675e-8im    -4.28056e-7-3.1921e-8im    …    4.01489e-7-1.44173e-8im
 -7.41517e-7-4.28364e-8im   -4.27914e-7-3.83248e-8im        4.01703e-7-1.7313e-8im  
 -7.41701e-7-5.00165e-8im   -4.27761e-7-4.47127e-8im        4.01912e-7-2.02371e-8im
 -7.41901e-7-5.72139e-8im   -4.27575e-7-5.11443e-8im        4.02202e-7-2.31354e-8im
 -7.42144e-7-6.44257e-8im   -4.27389e-7-5.75509e-8im        4.02492e-7-2.60687e-8im
 -7.42416e-7-7.16655e-8im   -4.27156e-7-6.39878e-8im   …    4.02823e-7-2.89802e-8im
 -7.42728e-7-7.89473e-8im   -4.26899e-7-7.04227e-8im        4.03202e-7-3.19224e-8im
 -7.43044e-7-8.62413e-8im   -4.26631e-7-7.68855e-8im        4.03602e-7-3.48864e-8im
            ⋮                                          ⋱                            
 -3.96623e-9-2.01634e-11im  -7.62991e-9-1.17542e-11im      8.40682e-10-2.82507e-12im
 -3.96938e-9+1.15232e-11im  -7.63317e-9+1.61364e-11im      8.21519e-10-1.31251e-11im
 -3.97244e-9-7.90479e-13im  -7.63068e-9+1.09033e-11im  …   8.11657e-10-5.99587e-12im
 -3.96541e-9+1.46034e-11im  -7.61813e-9-6.24478e-12im      8.26425e-10-6.05316e-12im
 -3.97574e-9+2.11076e-12im  -7.64141e-9+7.56728e-12im      8.29292e-10+1.46063e-11im
 -3.95737e-9-1.91491e-12im  -7.62251e-9+1.36544e-11im       8.3114e-10+1.63336e-12im
 -3.97089e-9+2.54508e-12im  -7.63469e-9-1.97353e-12im      8.57792e-10+8.34888e-14im
 -3.97389e-9+6.95222e-13im  -7.62701e-9-1.36602e-11im  …    8.2856e-10+1.14622e-11im
  -3.9615e-9+3.86602e-12im  -7.62552e-9+1.00606e-11im      8.41119e-10+1.14961e-11im
  -3.9738e-9+2.04636e-12im  -7.64419e-9+1.50613e-12im       8.5253e-10-3.44774e-12im
 -3.96815e-9-8.12039e-12im  -7.63589e-9+2.64089e-12im       8.4199e-10-4.30272e-12im
  -3.9762e-9+5.09676e-12im   -7.6402e-9-2.81547e-12im      8.31271e-10+3.22448e-12im
 -3.97375e-9+0.0im          -7.62793e-9+0.0im          …   8.40828e-10+0.0im      
```
!!! note

      By convention in Julia, data are stored *column-wise*. Here, `FFT` contains 188 individual Fourier transforms.

```@docs
FFTData
```
