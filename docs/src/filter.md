# `Filters` - methods for filtering `SeisData` and `SeisChannel` objects.

SeisNoise.jl provides `bandpass`, `bandstop`, `lowpass` and `highpass` filters built from
[DSP.jl](https://github.com/JuliaDSP/DSP.jl). Due to [multiple dispatch](https://en.wikipedia.org/wiki/Multiple_dispatch) in Julia, filter functions in SeisNoise.jl work on either the data in `SeisChannel` or `SeisData` objects or directly with Julia `Array`s. SeisNoise.jl uses Butterworth filters with a default of 4 corners.

## Filtering a SeisChannel or SeisData
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

## Filtering Arrays

Filtering arrays works much the same as filtering `SeisData` or `SeisChannel` objects.
The only additional variable required is the sampling rate `fs` of the data in the
array.

```julia
julia> cc_step, cc_len = 100, 100 # step and length of slices
julia> A, starts, ends = slide(S[1], cc_len, cc_step) # create sliding array
julia> fs = S[1].fs # get sampling rate
julia> demean!(A) # remove mean
julia> detrend!(A) # remove linear trend
julia> taper!(A,fs) # taper - defaults to 5% taper on either side of the trace
julia> bandpass(A,freqmin,freqmax,fs,corners=corners,zerophase=zerophase)
```




```@docs
bandpass!
bandstop!
lowpass!
highpass!
```
