# `Pre-Processing` - methods for cleaning raw noise data.

The pre-processing functions get raw data ready for correlation. Many of these
rely on [SeisIO's processing](https://seisio.readthedocs.io/en/latest/src/Processing/processing.html). A typical workflow includes the following steps:

- Loading data
- Filling gaps
- Downsampling
- Slicing day-long traces into smaller windows
- Demeaning, detrending and tapering

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

## Creating Sliding Windows
[Seats et al., 2011](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-246X.2011.05263.x) showed that short term stacks of cross-correlations converge more quickly when dividing raw data into short, overlapping windows. The `slide` function splits `S` into windows of length `cc_len` (in seconds), with step between windows `cc_step` (in seconds).

```julia
julia> cc_step, cc_len = 450, 1800 # 30 minute window length with 75% overlap
julia> A, starts, ends = slide(S[1], cc_len, cc_step)
julia> A
36000×189 Array{Float64,2}:
    0.0          -103.095  -67.9792   35.5011   -109.735   …   290.23     18.2172    2827.49    -8198.7      
    3.14758e-6   -105.382  -69.2759   29.9752   -105.452       483.704    80.7089   -2596.52    -9547.97     
    0.000105115  -104.929  -71.0974   24.1595   -104.596       635.558     9.51002  -1728.5     -7885.31     
    0.000551506  -104.545  -73.0085   16.088    -100.576      1193.93     30.3235   -4105.61    -3688.67     
    0.00155024   -104.646  -75.5247    6.76857   -95.9664     1181.32    -13.5895    1410.62      -65.8756   
    ⋮                                                      ⋱     ⋮                                           
 -104.322        -105.706  -54.9128  -65.9102     51.2037  …   733.128  -422.694      -76.919      -5.99003  
 -106.025        -106.238  -63.692   -44.0787     46.1554      155.186  -378.362      -63.3664     -1.61041  
 -108.454        -106.888  -70.0719  -26.3757     42.7611     -367.41   -310.284     -142.125       0.503033
 -113.258        -104.777  -79.5973  -20.8923     38.9552      180.431  -272.636     -137.636       0.167523
 -113.7          -103.16   -85.4747   -7.70531    36.072       472.996  -293.358      -74.2375     -0.0133164
```
Now `A` is a 2d Array containing 189 sliding windows, each 1800 seconds long.

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

## Pre-Processing Functions

```@docs
process_raw!
downsample
slide
detrend!
demean!
```
