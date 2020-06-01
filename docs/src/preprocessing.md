# Pre-Processing

The pre-processing functions get raw data ready for correlation. Many of these
rely on [SeisIO's processing](https://seisio.readthedocs.io/en/latest/src/Processing/processing.html). A typical workflow includes the following steps:

- Loading data
- Filling gaps
- Downsampling
- Slicing day-long traces into smaller windows
- Demeaning, detrending and tapering
- Time-domain normalization

Here is an example workflow:

## Load or download data

Data can be downloaded using `SeisIO`'s `get_data` function.

```julia
julia> using SeisIO, SeisNoise
julia> S = get_data("IRIS","TA.V04C..BHZ",s="2006-02-01",t="2006-02-02")
SeisData with 1 channels (1 shown)
    ID: TA.V04C..BHZ                       
  NAME: TA.V04C..BHZ                       
   LOC: 0.0 N, 0.0 E, 0.0 m                
    FS: 40.0                               
  GAIN: 1.0                                
  RESP: c = 0.0, 0 zeros, 0 poles          
 UNITS:                                    
   SRC:                                    
  MISC: 0 entries                          
 NOTES: 0 entries                          
     T: 2006-02-01T00:00:00.020 (0 gaps)   
     X: +1.044e-06                         
        +1.060e-06                         
            ...                            
        +6.282e-07                         
        (nx = 3456000)                     
     C: 0 open, 0 total
julia> writesac(S)
```

Data can be read from disk using the `read_data` function .

```julia
julia> S = read_data("sac","2006.032.00.00.00.002.TA.V04C..BHZ.R.SAC")
SeisData with 1 channels (1 shown)
    ID: TA.V04C..BHZ                       
  NAME:                                    
   LOC: 0.0 N, 0.0 E, 0.0 m                
    FS: 40.0                               
  GAIN: 1.0                                
  RESP: c = 0.0, 0 zeros, 0 poles          
 UNITS:                                    
   SRC: 2006.032.00.00.00.002.TA.V04C..BH…
  MISC: 0 entries                          
 NOTES: 1 entries                          
     T: 2006-02-01T00:00:00.002 (0 gaps)   
     X: +1.044e-06                         
        +1.060e-06                         
            ...                            
        +6.282e-07                         
        (nx = 3456000)                     
     C: 0 open, 0 total
```



## Merge and ungap data
Use the `merge!` method to merge `S` to have one channel. The `ungap` method replaces
gaps in `S[1].x` with the mean of `S[1].x`.

```julia
julia> merge!(S)
julia> ungap!(S)
```

## Downsample Data  
Downsample `S` to 20 Hz using `process_raw!`. `process_raw!`:
- Removes mean from each channel in `S`.
- Detrends each channel in `S`.
- Downsamples data to sampling rate `fs`
- Phase-shifts data to begin at 00:00:00.0
    - This is important for data that does not begin exactly on the sampling rate, e.g. the starttime is 2006-02-01T00:00:00.020.  

```julia
julia> process_raw!(S,20.)
julia> S
SeisData with 1 channels (1 shown)
    ID: TA.V04C..BHZ                       
  NAME: TA.V04C..BHZ                       
   LOC: 0.0 N, 0.0 E, 0.0 m                
    FS: 20.0                               
  GAIN: 1.0                                
  RESP: c = 0.0, 0 zeros, 0 poles          
 UNITS:                                    
   SRC:                                    
  MISC: 0 entries                          
 NOTES: 1 entries                          
     T: 2006-02-01T00:00:00.000 (0 gaps)   
     X: +9.289e-10                         
        -9.248e-10                         
            ...                            
        -9.287e-10                         
        (nx = 1728000)                     
     C: 0 open, 0 total
```

Note that now `S` has sampling rate `fs = 20.0`, has half as many points (`nx = 1728001`) as before and the starttime has changed to 2006-02-01T00:00:00.000.

## Detrending and Demeaning

The `demean` and `detrend` functions are applied column wise.

```julia
julia> A = reshape(collect(1.:10.),5,2)
5×2 Array{Float64,2}:
 1.0   6.0
 2.0   7.0
 3.0   8.0
 4.0   9.0
 5.0  10.0
julia> demean(A)
5×2 Array{Float64,2}:
 -2.0  -2.0
 -1.0  -1.0
  0.0   0.0
  1.0   1.0
  2.0   2.0
julia> detrend(A)
5×2 Array{Float64,2}:
6.66134e-16  3.55271e-15
1.33227e-15  5.32907e-15
2.22045e-15  7.10543e-15
2.66454e-15  7.10543e-15
3.55271e-15  1.06581e-14
```

`demean` and `detrend` also work on `RawData` and `CorrData` structure. If `R` is a `RawData`, then the windowed data stored in `R.x` can be detrended in-place using

```julia
detrend!(R)
```

or allocated to a new `RawData` using

```julia
Rd = detrend(R)
```

## Amplitude Normaliztion

Time-domain normalization is used to suppress high-amplitude signals, such as earthquakes or instrumental irregularities. SeisNoise provides time-normalization functions for:

- one-bit normalization: `onebit`
- clipping: `clip`
- running mean normalization: `smooth`
- suppressing high-amplitude signals: `mute`


## Filtering
SeisNoise.jl provides `bandpass`, `bandstop`, `lowpass` and `highpass` filters built from
[DSP.jl](https://github.com/JuliaDSP/DSP.jl). Due to [multiple dispatch](https://en.wikipedia.org/wiki/Multiple_dispatch) in Julia, filter functions in SeisNoise.jl work on either the data in `RawData` or `CorrData` objects or directly with Julia `Array`s. SeisNoise.jl uses Butterworth filters with a default of 4 corners.

### Filtering a RawData or CorrData
```julia
julia> freqmin, freqmax = 1., 10. # low and high frequencies in Hz
julia> corners = 4 # number of corners in Butterworth filter
julia> zerophase = true # if true, apply filter twice for no phase shifting
julia> S = get_data("IRIS","TA.V04C..BHZ",s="2006-02-01T00:00:00",t="2006-02-01T01:00:00")
julia> SeisIO.demean!(S) # remove mean
julia> SeisIO.detrend!(S) # remove linear trend
julia> SeisIO.taper!(S) # taper - defaults to 5% taper on either side of the trace
julia> bandpass(S,freqmin,freqmax,corners=corners,zerophase=zerophase)
```

### Nyquist Frequency
Filtering above the Nyquist frequency will give a warning, if using a `lowpass`
or `bandpass` filter, while using a `highpass` above the Nyquist frequency will throw an
error.

```julia
julia> S[1].fs / 2. # Nyquist frequency is half sampling rate
20.
julia> freqmax = S[1].fs / 2. + 5.
25.  
julia> bandpass(S,freqmin,freqmax,corners=corners,zerophase=zerophase)
┌ Warning: Selected high corner frequency (25.0) of bandpass is at or
│        above Nyquist (20.0). Applying a high-pass instead.
└
julia> lowpass(S,freqmax,corners=corners,zerophase=zerophase)
┌ Warning: Selected corner frequency (25.0) is
│ above Nyquist (20.0). Setting Nyquist as high corner.
└
julia> highpass(S,freqmax,corners=corners,zerophase=zerophase)
ERROR: frequencies must be less than the Nyquist frequency 20.0
```


## Pre-Processing Functions

```@autodocs
Modules = [SeisNoise]
Pages   = ["ArrayFuncs.jl", "filter.jl"]
```

```@docs
process_raw!
onebit!
clip!
clamp!
mute!
```
