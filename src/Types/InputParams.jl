export InputParams

"""
    InputParams

A structure for inputing parameters necessary for processing ambient noise data.

## Fields: InputParams
| **Field** | **Description** |
|:--------------------|:--------------------|
| :pair               | Network.Station.Channel for each station to be correlated. |
| :FFTOUT             | Directory filepath for FFT output. |
| :CORROUT            | Directory filepath for correlation output. |
| :cc_len             | Length of each correlation in seconds. |
| :cc_step            | Spacing between correlation windows in seconds. |
| :fs                 | Sampling frequency in Hz to downsample raw data; set to 0.0. |
| :freqmin            | Minimum frequency for whitening.  |
| :freqmax            | Maximum frequency for whitening. |
| :time_norm          | Apply one-bit whitening with "one_bit". |
| :to_whiten          | Apply whitening to FFTs. |
| :maxlag             | Maximum lag time in seconds to keep from correlations. |
| :smoothing_half_win | Half window in samples for FFT smoothing. |
| :corr_type          | Type of correlation e.g. "cross-correlate", "deconv" or "coherence" |
"""
mutable struct InputParams
           pair::Array{String,1}
           FFTOUT::String
           CORROUT::String
           cc_len::Int
           cc_step::Int
           fs::Float64
           freqmin::Float64
           freqmax::Float64
           time_norm::Union{Bool,String}
           to_whiten::Bool
           maxlag::Float64
           smoothing_half_win::Int
           corr_type::String

     function InputParams(
           pair::Array{String,1},
           FFTOUT::String,
           CORROUT::String,
           cc_len::Int,
           cc_step::Int,
           fs::Float64,
           freqmin::Float64,
           freqmax::Float64,
           time_norm::Union{Bool,String},
           to_whiten::Bool,
           maxlag::Float64,
           smoothing_half_win::Int,
           corr_type::String)

           return new(pair,FFTOUT,CORROUT,cc_len,cc_step,fs,freqmin,freqmax,
                      time_norm, to_whiten, maxlag, smoothing_half_win,
                      corr_type)
           end
end

InputParams(;
            pair::Array{String,1}         = Array{String,1}(undef, 0),
            FFTOUT::String                = "",
            CORROUT::String               = "",
            cc_len::Int                   = zero(Int),
            cc_step::Int                  = zero(Int),
            fs::Float64                   = zero(Float64),
            freqmin::Float64              = zero(Float64),
            freqmax::Float64              = zero(Float64),
            time_norm::Union{Bool,String} = false,
            to_whiten::Bool               = false,
            maxlag::Float64               = zero(Float64),
            smoothing_half_win::Int       = zero(Int),
            corr_type::String             = ""
            ) = InputParams(pair,FFTOUT,CORROUT,cc_len,cc_step,fs,freqmin,freqmax,
                       time_norm, to_whiten, maxlag, smoothing_half_win,
                       corr_type)
