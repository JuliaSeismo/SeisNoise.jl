var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Noise.jl-1",
    "page": "Home",
    "title": "Noise.jl",
    "category": "section",
    "text": "Ambient Noise Cross-Correlation in Julia.Noise.jl provides routines for quickly and efficiently implementing seismic interferometry.note: Note\nMuch of the I/O and preprocessing in Noise.jl uses the SeisIO.jl package. Please read through the SeisIO Documentation to get familiar with seismic data processing in Julia."
},

{
    "location": "#Package-Features-1",
    "page": "Home",
    "title": "Package Features",
    "category": "section",
    "text": "Built upon SeisIO for easy and fast I/O.\nCustom types for saving Fourier Transforms of data and cross-correlations\nArray-based processing of raw data and cross-correlation.\nComing soon: methods for dvv and dispersion measurements.\nComing soon: GPU support."
},

{
    "location": "#Processing-1",
    "page": "Home",
    "title": "Processing",
    "category": "section",
    "text": ""
},

{
    "location": "#Plotting-1",
    "page": "Home",
    "title": "Plotting",
    "category": "section",
    "text": ""
},

{
    "location": "#High-Performance-Computing-1",
    "page": "Home",
    "title": "High-Performance Computing",
    "category": "section",
    "text": ""
},

{
    "location": "#Tutorial-1",
    "page": "Home",
    "title": "Tutorial",
    "category": "section",
    "text": ""
},

{
    "location": "preprocessing/#",
    "page": "Pre-Processing",
    "title": "Pre-Processing",
    "category": "page",
    "text": ""
},

{
    "location": "preprocessing/#Noise.process_raw!",
    "page": "Pre-Processing",
    "title": "Noise.process_raw!",
    "category": "function",
    "text": "process_raw!(S,fs)\n\nPre-process raw seismic data.\n\nChecks:\n\nsample rate is fs\ndownsamples data\nchecks for gaps in data\nphase-shifts data to begin at 00:00:00.0\n\nArguments\n\nS::SeisChannel: SeisData structure.\nfs::Float64: Sampling rate to downsample S.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Noise.downsample",
    "page": "Pre-Processing",
    "title": "Noise.downsample",
    "category": "function",
    "text": "downsample(C::SeisChannel,fs::Real)\n\nDownsample SeisChannel sampling rate to frequency fs.\n\nFor best results, lowpass filter data to fs before downsampling. Implements the weighted average slopes interpolation scheme proposed in [Wiggins1976] for evenly sampled data from obspy.signal.interpolation.weightedaverageslopes.\n\n\n\n\n\ndownsample(S::SeisData, fs::Float64)\n\nDownsample SeisData sampling rate to frequency fs.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Noise.slide",
    "page": "Pre-Processing",
    "title": "Noise.slide",
    "category": "function",
    "text": "slide(C::SeisChannel, window_length::Real, cc_step::Real)\n\nGenerate equal length sliding windows into an array.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Noise.start_end",
    "page": "Pre-Processing",
    "title": "Noise.start_end",
    "category": "function",
    "text": "start_end(C::SeisChannel)\n\nReturn start and endtimes of SeisChannel in DateTime format.\n\n\n\n\n\nstart_end(S::SeisData)\n\nReturn start and endtimes of SeisChannel in DateTime format.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Noise.nearest_start_end",
    "page": "Pre-Processing",
    "title": "Noise.nearest_start_end",
    "category": "function",
    "text": "nearest_start_end(C::SeisChannel, cc_len::Int, cc_step::Int)\n\nReturn best possible start, end times for data in C given the c\n\n\n\n\n\nnearest_start_end(D::DateTime, cc_len::Int, cc_step::Int)\n\nReturn best possible start, end times for given starttime D\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Pre-Processing-methods-for-cleaning-raw-noise-data.-1",
    "page": "Pre-Processing",
    "title": "Pre-Processing - methods for cleaning raw noise data.",
    "category": "section",
    "text": "process_raw!\ndownsample\nslide\nstart_end\nnearest_start_end"
},

{
    "location": "arrayfuncs/#",
    "page": "Array Functions",
    "title": "Array Functions",
    "category": "page",
    "text": ""
},

{
    "location": "arrayfuncs/#Noise.detrend!",
    "page": "Array Functions",
    "title": "Noise.detrend!",
    "category": "function",
    "text": "detrend!(X::AbstractArray{<:Union{Float32,Float64},1})\n\nRemove linear trend from array X using least-squares regression.\n\n\n\n\n\ndetrend!(X::AbstractArray{<:Union{Float32,Float64},2})\n\nRemove linear trend from columns of X using least-squares regression.\n\n\n\n\n\n"
},

{
    "location": "arrayfuncs/#Noise.demean!",
    "page": "Array Functions",
    "title": "Noise.demean!",
    "category": "function",
    "text": "demean!(A::AbstractArray{<:Union{Float32,Float64},1})\n\nRemove mean from array A.\n\n\n\n\n\ndemean!(A::AbstractArray{<:Union{Float32,Float64},2})\n\nRemove mean from columns of array A.\n\n\n\n\n\n"
},

{
    "location": "arrayfuncs/#Array-Functions-methods-for-operating-on-array-based-noise-and-corelations.-1",
    "page": "Array Functions",
    "title": "Array Functions - methods for operating on array-based noise and corelations.",
    "category": "section",
    "text": "detrend!\ndemean!"
},

{
    "location": "filter/#",
    "page": "Filters",
    "title": "Filters",
    "category": "page",
    "text": ""
},

{
    "location": "filter/#Noise.bandpass!",
    "page": "Filters",
    "title": "Noise.bandpass!",
    "category": "function",
    "text": "bandpass!(A,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandpass Filter.\n\nFilter data A from freqmin to freqmax using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreqmin::Float64: Pass band low corner frequency.\nfreqmax::Float64: Pass band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nbandpass!(C,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandpass Filter.\n\nFilter data in C from freqmin to freqmax using corners corners.\n\nArguments\n\nC::SeisChannel: SeisChannel to filter.\nfreqmin::Float64: Pass band low corner frequency.\nfreqmax::Float64: Pass band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nbandpass!(S,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandpass Filter.\n\nFilter channels in S from freqmin to freqmax using corners corners.\n\nArguments\n\nS::SeisData: SeisData to filter.\nfreqmin::Float64: Pass band low corner frequency.\nfreqmax::Float64: Pass band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "filter/#Noise.bandstop!",
    "page": "Filters",
    "title": "Noise.bandstop!",
    "category": "function",
    "text": "bandstop!(A,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandstop Filter.\n\nFilter data A removing data between frequencies freqmin to freqmax using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreqmin::Float64: Stop band low corner frequency.\nfreqmax::Float64: Stop band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nbandstop!(C,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandstop Filter.\n\nFilter data in C removing data between frequencies freqmin to freqmax using corners corners.\n\nArguments\n\nC::SeisChannel: SeisChannel to filter.\nfreqmin::Float64: Stop band low corner frequency.\nfreqmax::Float64: Stop band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "filter/#Noise.lowpass!",
    "page": "Filters",
    "title": "Noise.lowpass!",
    "category": "function",
    "text": "lowpass(A,freq,fs,corners=4,zerophase=false)\n\nButterworth-Lowpass Filter.\n\nFilter data A over certain frequency freq using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nlowpass(C,freq,fs,corners=4,zerophase=false)\n\nButterworth-Lowpass Filter.\n\nFilter data in C over certain frequency freq using corners corners.\n\nArguments\n\nC::SeisChannel: SeisChannel to filter.\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nlowpass(S,freq,fs,corners=4,zerophase=false)\n\nButterworth-Lowpass Filter.\n\nFilter channels in S over certain frequency freq using corners corners.\n\nArguments\n\nS::SeisData: SeisData to filter.\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "filter/#Noise.highpass!",
    "page": "Filters",
    "title": "Noise.highpass!",
    "category": "function",
    "text": "highpass(A,freq,fs,corners=4,zerophase=false)\n\nButterworth-Highpass Filter.\n\nFilter data A removing data below certain frequency freq using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nhighpass(C,freq,fs,corners=4,zerophase=false)\n\nButterworth-Highpass Filter.\n\nFilter data in C removing data below certain frequency freq using corners corners.\n\nArguments\n\nC::SeisChannel: SeisChannel to filter.\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nhighpass(S,freq,fs,corners=4,zerophase=false)\n\nButterworth-Highpass Filter.\n\nFilter channels in S removing data below certain frequency freq using corners corners.\n\nArguments\n\nS::SeisData: SeisData to filter.\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "filter/#Filters-methods-for-filtering-SeisData-and-SeisChannel-objects.-1",
    "page": "Filters",
    "title": "Filters - methods for filtering SeisData and SeisChannel objects.",
    "category": "section",
    "text": "bandpass!\nbandstop!\nlowpass!\nhighpass!"
},

{
    "location": "fft/#",
    "page": "Computing FFTs",
    "title": "Computing FFTs",
    "category": "page",
    "text": ""
},

{
    "location": "fft/#Noise.whiten",
    "page": "Computing FFTs",
    "title": "Noise.whiten",
    "category": "function",
    "text": "whiten(A, freqmin, freqmax, fs, pad=100)\n\nWhiten spectrum of time series A between frequencies freqmin and freqmax. Uses real fft to speed up computation. Returns the whitened (single-sided) fft of the time series.\n\nArguments\n\nA::AbstractArray: Time series.\nfs::Real: Sampling rate of time series A.\nfreqmin::Real: Pass band low corner frequency.\nfreqmax::Real: Pass band high corner frequency.\npad::Int: Number of tapering points outside whitening band.\n\n\n\n\n\n"
},

{
    "location": "fft/#Noise.process_fft",
    "page": "Computing FFTs",
    "title": "Noise.process_fft",
    "category": "function",
    "text": "process_fft(A::AbstractArray,freqmin::Float64,freqmax::Float64,fs::Float64;\n            time_norm=false,to_whiten=false,corners=corners,\n            zerophase=zerophase)\n\nArguments\n\nA::AbstractArray: Array with time domain data.\nfs::Float64: Sampling rate of data in A.\nfreqmin::Float64: minimum frequency for whitening.\nfreqmax::Float64: maximum frequency for whitening.\ntime_norm::Union{Bool,String}: time domain normalization to perform.\nto_whiten::Bool: Apply whitening in frequency domain.\ncorners::Int: Number of corners in Butterworth filter.\nzerophase::Bool: If true, apply Butterworth filter twice for zero phase                    change in output signal.\n\n\n\n\n\n"
},

{
    "location": "fft/#Noise.save_fft",
    "page": "Computing FFTs",
    "title": "Noise.save_fft",
    "category": "function",
    "text": "save_fft(F::FFTData, OUT::String)\n\nSave FFTData F to JLD2.\n\n\n\n\n\n"
},

{
    "location": "fft/#Computing-FFTs-methods-for-computing-FFTs-raw-noise-data.-1",
    "page": "Computing FFTs",
    "title": "Computing FFTs - methods for computing FFTs raw noise data.",
    "category": "section",
    "text": "whiten\nprocess_fft\nsave_fft"
},

{
    "location": "correlation/#",
    "page": "Correlation",
    "title": "Correlation",
    "category": "page",
    "text": ""
},

{
    "location": "correlation/#Noise.clean_up!",
    "page": "Correlation",
    "title": "Noise.clean_up!",
    "category": "function",
    "text": "clean_up!(A,freqmin,freqmax,fs)\n\nDemean, detrend, taper and filter time series.\n\nArguments\n\nA::AbstractArray: Time series.\nfs::Real: Sampling rate of time series A.\nfreqmin::Real: Pass band low corner frequency.\nfreqmax::Real: Pass band high corner frequency.\n\n\n\n\n\n"
},

{
    "location": "correlation/#Noise.correlate",
    "page": "Correlation",
    "title": "Noise.correlate",
    "category": "function",
    "text": "correlate(fft1, fft2, N, maxlag, corr_type=\'cross-correlate\')\n\nCross-correlation of two ffts.\n\n\n\n\n\n"
},

{
    "location": "correlation/#Noise.compute_cc",
    "page": "Correlation",
    "title": "Noise.compute_cc",
    "category": "function",
    "text": "compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;\n           smoothing_half_win::Int=20,\n           corr_type::String=\"cross-correlation\" )\n\n\n\n\n\n"
},

{
    "location": "correlation/#Noise.next_fast_len",
    "page": "Correlation",
    "title": "Noise.next_fast_len",
    "category": "function",
    "text": "next_fast_len(N::Real)\n\nReturn next fast length for fft with FFTW.\n\n\n\n\n\n"
},

{
    "location": "correlation/#Noise.save_corr",
    "page": "Correlation",
    "title": "Noise.save_corr",
    "category": "function",
    "text": "save_corr(C::CorrData, OUT::String)\n\nSave CorrData C to JLD2.\n\n\n\n\n\n"
},

{
    "location": "correlation/#Noise.load_fft",
    "page": "Correlation",
    "title": "Noise.load_fft",
    "category": "function",
    "text": "load_fft(filename,chan,day)\n\nLoads FFTData for channel chan on day day from JLD2 file filename.\n\n\n\n\n\n"
},

{
    "location": "correlation/#Computing-Correlations-methods-for-computing-correlations-from-FFTs.-1",
    "page": "Correlation",
    "title": "Computing Correlations - methods for computing correlations from FFTs.",
    "category": "section",
    "text": "clean_up!\ncorrelate\ncompute_cc\nnext_fast_len\nsave_corr\nload_fft"
},

{
    "location": "Types/inputparams/#Noise.InputParams",
    "page": "InputParams",
    "title": "Noise.InputParams",
    "category": "type",
    "text": "InputParams\n\nA structure for inputing parameters necessary for processing ambient noise data.\n\nFields: InputParams\n\nField Description\n:pair Network.Station.Channel for each station to be correlated.\n:FFTOUT Directory filepath for FFT output.\n:CORROUT Directory filepath for correlation output.\n:cc_len Length of each correlation in seconds.\n:cc_step Spacing between correlation windows in seconds.\n:fs Sampling frequency in Hz to downsample raw data; set to 0.0.\n:freqmin Minimum frequency for whitening.\n:freqmax Maximum frequency for whitening.\n:time_norm Apply one-bit whitening with \"one_bit\".\n:to_whiten Apply whitening to FFTs.\n:maxlag Maximum lag time in seconds to keep from correlations.\n:smoothinghalfwin Half window in samples for FFT smoothing.\n:corr_type Type of correlation e.g. \"cross-correlate\", \"deconv\" or \"coherence\"\n\n\n\n\n\n"
},

{
    "location": "Types/inputparams/#",
    "page": "InputParams",
    "title": "InputParams",
    "category": "page",
    "text": "Input Parameters - Objects for scheduling cross-correlationsInputParams"
},

{
    "location": "Types/fftdata/#Noise.FFTData",
    "page": "FFTData",
    "title": "Noise.FFTData",
    "category": "type",
    "text": "FFTData\n\nA structure for fourier transforms (FFT) of ambient noise data.\n\nFields: FFTData\n\nField Description\n:name Freeform channel names\n:id Channel ids. use NET.STA.LOC.CHAN format when possible.\n:loc Location (position) vector; freeform.\n:fs Sampling frequency in Hz.\n:gain Scalar gain; divide data by the gain to convert to units\n:freqmin Minimum frequency for whitening.\n:freqmax Maximum frequency for whitening.\n:cc_len Length of each correlation in seconds.\n:cc_step Spacing between correlation windows in seconds.\n:whitened Whitening applied.\n:time_norm Apply one-bit whitening with \"one_bit\".\n:resp Instrument response; two-column matrix, format [zeros poles]\n:misc Dictionary for non-critical information.\n:notes Timestamped notes; includes automatically-logged acquisition and\n processing information.\n:maxlag Maximum lag time in seconds to keep from correlations.\n:t Starttime of each FFT.\n:fft FFTs stored in columns.\n\n\n\n\n\n"
},

{
    "location": "Types/fftdata/#",
    "page": "FFTData",
    "title": "FFTData",
    "category": "page",
    "text": "FFTData - Objects for holding Fourier transforms (FFTs)FFTData"
},

{
    "location": "Types/corrdata/#Noise.CorrData",
    "page": "CorrData",
    "title": "Noise.CorrData",
    "category": "type",
    "text": "CorrData\n\nA structure for cross-correlations of ambient noise data.\n\nFields: CorrData\n\nField Description\n:name Freeform channel names\n:id Channel ids. use NET.STA.LOC.CHAN format when possible.\n:loc Location (position) vector; freeform.\n:fs Sampling frequency in Hz.\n:gain Scalar gain; divide data by the gain to convert to units\n:freqmin Minimum frequency for whitening.\n:freqmax Maximum frequency for whitening.\n:cc_len Length of each correlation in seconds.\n:cc_step Spacing between correlation windows in seconds.\n:whitened Whitening applied.\n:time_norm Apply one-bit whitening with \"one_bit\".\n:resp Instrument response; two-column matrix, format [zeros poles]\n:misc Dictionary for non-critical information.\n:notes Timestamped notes; includes automatically-logged acquisition and\n processing information.\n:maxlag Maximum lag time in seconds to keep from correlations.\n:t Starttime of each correlation.\n:corr Correlations stored in columns.\n\n\n\n\n\n"
},

{
    "location": "Types/corrdata/#",
    "page": "CorrData",
    "title": "CorrData",
    "category": "page",
    "text": "CorrData - Objects for ambient noise cross-correlations.CorrData"
},

]}
