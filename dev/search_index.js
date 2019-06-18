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
    "location": "correlation/#Computing-Correlations-methods-for-computing-correlations-from-FFTs.-1",
    "page": "Correlation",
    "title": "Computing Correlations - methods for computing correlations from FFTs.",
    "category": "section",
    "text": "clean_up!\ncorrelate\ncompute_cc\nnext_fast_len\nsave_corr\nload_fft"
},

{
    "location": "postprocessing/#",
    "page": "Velocity Change",
    "title": "Velocity Change",
    "category": "page",
    "text": ""
},

{
    "location": "postprocessing/#Noise.stack!",
    "page": "Velocity Change",
    "title": "Noise.stack!",
    "category": "function",
    "text": "stack!(C)\n\nStack correlation by time interval. The default is to stack by day. Using allstack == true will stack all available correlations. To use phase-weighted stack, specify the amount of phase_smoothing in seconds.\n\nArguments\n\nC::CorrData: Correlation data.\ninterval::Union{Month,Day,Hour,Second}: Interval over which to stack C.\nallstack::Bool: If true, stack all data.\nphase_smoothing::Float64: Enables phase-weighted stacking. phase_smoothing                             is the time window in seconds for phase smoothing                             in the phase-weighted stack.\n\n\n\n\n\n"
},

{
    "location": "postprocessing/#Noise.mwcs",
    "page": "Velocity Change",
    "title": "Noise.mwcs",
    "category": "function",
    "text": "mwcs(ref, cur, fmin, fmax, fs, tmin, windowlength, windowstep,        smoothinghalfwin)\n\nChange in velocity measurement using the Moving Window Cross-Spectrum technique.\n\nThe current correlation cur is compared to the reference correlation ref. Both time series are sliced in several overlapping windows. Each slice is mean-adjusted and cosine-tapered (85% taper) before being Fourier- transformed to the frequency domain. F_cur(ν) and F_ref(ν) are the first halves of the Hermitian symmetric Fourier-transformed segments. The cross-spectrum `X(ν) is defined as X(ν) = F_ref(ν) F_cur^*(ν) in which ^* denotes the complex conjugation. X(ν) is then smoothed by convolution with a Hanning window. The similarity of the two time-series is assessed using the cross-coherency between energy densities in the frequency domain: C(ν) = fracoverlineX(ν))sqrtoverlineF_ref(ν)^2 overlineF_cur(ν)^2 in which the over-line here represents the smoothing of the energy spectra for F_ref and F_cur and of the spectrum of X. The mean coherence for the segment is defined as the mean of C(ν) in the frequency range of interest. The time-delay between the two cross correlations is found in the unwrapped phase, ϕ(ν), of the cross spectrum and is linearly proportional to frequency: ϕ_j = m ν_j m = 2 π δ t The time shift for each window between two signals is the slope m of a weighted linear regression of the samples within the frequency band of interest. The weights are those introduced by [Clarke2011], which incorporate both the cross-spectral amplitude and cross-coherence, unlike [Poupinet1984]. The errors are estimated using the weights (thus the coherence) and the squared misfit to the modelled slope: e_m = sqrtsum_j(fracw_j ν_jsum_iw_i ν_i^2)^2σ_ϕ^2 where w are weights, ν are cross-coherences and σ_ϕ^2 is the squared misfit of the data to the modelled slope and is calculated as σ_ϕ^2 = fracsum_j(ϕ_j - m ν_j)^2N-1 The output of this process is a table containing, for each moving window: the central time lag, the measured delay, its error and the mean coherence of the segment.\n\nArguments\n\nref::AbstractArray: Reference correlation.\ncur::AbstractArray: Current correlation.\nfmin:Float64: minimum frequency in the correlation [Hz]\nfmax:Float64: maximum frequency in the correlation [Hz]\nfs:Float64: Sampling rate of ref and cur [Hz]\ntmin:Float64: The leftmost time lag [s]\nwindow_len:Float64: The moving window length [s]\nwindow_step:Float64: The step to jump for the moving window [s]\nsmoothing_half_win:Int: Defines the half length of the smoothing hanning                           window.\n\nReturns\n\ntime_axis:Array{Float64,1}: Central time of each moving window [s]\ndt:Array{Float64,1}: dt for each moving window\nerr:Array{Float64,1}: Errors for each moving window\nmcoh:Array{Float64,1}: Mean coherence for each moving window\n\nThis code is a Julia translation of the Python code from MSNoise.\n\n\n\n\n\n"
},

{
    "location": "postprocessing/#Noise.mwcs_dvv",
    "page": "Velocity Change",
    "title": "Noise.mwcs_dvv",
    "category": "function",
    "text": "mwcs_dvv(time_axis,dt,err,coh,dtt_lag,dist,dtt_v,dtt_minlag,dtt_width,\n         dtt_sides,min_coh, max_err, max_dt)\n\nRegresses dv/v from dt/t measurements.\n\nSolves dt = a + bt, where b = dv/v, a = instrumental drift, dt = time lag at time t. Solves with a weighted linear regression with weights equal to 1/error**2.\n\nArguments\n\ntime_axis:Array{Float64,1}: Central time of each moving window [s]\ndt:Array{Float64,1}: dt for each moving window\nerr:Array{Float64,1}: Errors for each moving window\ncoh:Array{Float64,1}: Mean coherence for each moving window\ndist:Float64: Distance between stations [km]\ndtt_lag:String: Type of time lag \'dynamic\' or \'static\'. When dtt_lag   is set to \"dynamic\", the inter-station distance is used to determine   the minimum time lag\ndtt_v:Float64: Velocity for minumum time lag. The velocity is determined by   the user so that the minlag doesn\'t include the ballistic waves.   If ballistic waves are visible with a velocity of 2 km/s, one could   configure dtt_v=1.0. If stations are located 15 km apart, the minimum   lag time will be set to 15 s.\nminlag:Float64: Statically set minimum lag [s]\ndtt_width:Float64: Width of the lag window used [s] i.e., if   dttwidth = 30, and minlag = 15, window = 15s - 45s\ndtt_sides:String: Sides of cross-correlation function to use. Either   \'Both\' or \'left\'.\nmin_coh:Float64: Minimum allowed coherency between reference and current                    correlation\nmax_err:Float64: Maximum allowed error in dt/t regression\nmax_dt:Float64: Maximum allowed dt from MWCS [s]\n\nReturns\n\nm::Float64: dt/t for current correlation\nem::Float64: Error for calulatoin of m\na::Float64: Intercept for regression calculation\nea::Float64: Error on intercept\nm0::Float64: dt/t for current correlation with no intercept\nem0::Float64: Error for calulatoin of m0\n\n\n\n\n\n"
},

{
    "location": "postprocessing/#Noise.stretching",
    "page": "Velocity Change",
    "title": "Noise.stretching",
    "category": "function",
    "text": "stretching(ref,cur,t,window,fmin,fmax;dvmin,dvmax,ntrials)\n\nThis function compares the Reference waveform to stretched/compressed current waveforms to get the relative seismic velocity variation (and associated error). It also computes the correlation coefficient between the Reference waveform and the current waveform.\n\nArguments\n\nref::AbstractArray: Reference correlation.\ncur::AbstractArray: Current correlation.\nt::AbstractArray: time vector, common to both ref and cur.\nwindow:AbstractArray: vector of the indices of the cur and ref windows                         on which you want to do the measurements\nfmin:Float64: minimum frequency in the correlation [Hz]\nfmax:Float64: maximum frequency in the correlation [Hz]\ndvmin:Float64: minimum bound for the velocity variation; e.g. dvmin=-0.03                  for -3% of relative velocity change\ndvmax:Float64: maximum bound for the velocity variation; e.g. dvmin=0.03                 for 3% of relative velocity change\nntrial:  number of stretching coefficient between dvmin and dvmax, no need to be higher than 100\n\nReturns\n\ndv:AFloat64: Relative Velocity Change dv/v (in %)\ncc:Float64: Correlation coefficient between the reference waveform and the                     best stretched/compressed current waveform\ncdp:Float64: Correlation coefficient between the reference waveform and the                initial current waveform\nϵ:Array{Float64,1}: Vector of Epsilon values (ϵ =-dt/t = dv/v)\nerr:Float64: Errors in the dv/v measurements based on Weaver et al., 2011\nC:Array{Float64,1}: Vector of the correlation coefficient between the                       reference waveform and every stretched/compressed                       current waveforms\n\nThis code is a Julia translation of the Python code from Viens et al., 2018.\n\n\n\n\n\n"
},

{
    "location": "postprocessing/#Post-Processing-methods-for-working-with-correlation-data.-1",
    "page": "Velocity Change",
    "title": "Post-Processing - methods for working with correlation data.",
    "category": "section",
    "text": "stack!\nmwcs\nmwcs_dvv\nstretching"
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
    "text": "FFTData\n\nA structure for fourier transforms (FFT) of ambient noise data.\n\nFields: FFTData\n\nField Description\n:name Freeform channel names\n:id Channel ids. use NET.STA.LOC.CHAN format when possible.\n:loc Location (position) object\n:fs Sampling frequency in Hz.\n:gain Scalar gain; divide data by the gain to convert to units\n:freqmin Minimum frequency for whitening.\n:freqmax Maximum frequency for whitening.\n:cc_len Length of each correlation in seconds.\n:cc_step Spacing between correlation windows in seconds.\n:whitened Whitening applied.\n:time_norm Apply one-bit whitening with \"one_bit\".\n:resp Instrument response object, format [zeros poles]\n:misc Dictionary for non-critical information.\n:notes Timestamped notes; includes automatically-logged acquisition and\n processing information.\n:t Starttime of each FFT.\n:fft FFTs stored in columns.\n\n\n\n\n\n"
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
    "text": "CorrData\n\nA structure for cross-correlations of ambient noise data.\n\nFields: CorrData\n\nField Description\n:name Freeform channel names\n:id Channel ids. use NET.STA.LOC.CHAN format when possible.\n:loc Location (position) object.\n:fs Sampling frequency in Hz.\n:gain Scalar gain; divide data by the gain to convert to units\n:freqmin Minimum frequency for whitening.\n:freqmax Maximum frequency for whitening.\n:cc_len Length of each correlation in seconds.\n:cc_step Spacing between correlation windows in seconds.\n:whitened Whitening applied.\n:time_norm Apply one-bit whitening with \"one_bit\".\n:resp Instrument response object, format [zeros poles]\n:misc Dictionary for non-critical information.\n:notes Timestamped notes; includes automatically-logged acquisition and\n processing information.\n:maxlag Maximum lag time in seconds to keep from correlations.\n:t Starttime of each correlation.\n:corr Correlations stored in columns.\n\n\n\n\n\n"
},

{
    "location": "Types/corrdata/#",
    "page": "CorrData",
    "title": "CorrData",
    "category": "page",
    "text": "CorrData - Objects for ambient noise cross-correlations.CorrData"
},

]}
