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
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "From the Julia command prompt:Press ] to enter pkg.\nType or copy: add https://github.com/tclements/Noise.jl; build; precompile\nPress backspace to exit pkg.\nType or copy: using Noise"
},

{
    "location": "#Using-Noise.jl-1",
    "page": "Home",
    "title": "Using Noise.jl",
    "category": "section",
    "text": "Noise.jl was designed to be as easy to use in the REPL as on an HPC cluster. If you want to get started processing data, head over to the tutorial or parallel example. This documentation provides a reference to all the underlying function for cross-correlation. We encourage you to develop your own workflow using Noise\'s core functionality.(Image: plot1)"
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
    "text": "FFTData - Objects for holding Fourier transforms (FFTs)Noise is designed around array-based cross-correlation. Noise uses a custom structure FFTData for holding Fourier transforms of ambient noise. To create an empty FFTData object, use the FFTData() function. FFTData are created with the compute_fft function.julia> using Noise\njulia> FFTData()\nFFTData with 0 ffts\n      NAME: …\n        ID: …\n       LOC: 0.0 N, 0.0 E, 0.0 m\n        FS: 0.0\n      GAIN: 1.0\n   FREQMIN: 0.0\n   FREQMAX: 0.0\n    CC_LEN: …\n   CC_STEP: …\n  WHITENED: …\n TIME_NORM: …\n      RESP: c = 0.0, 0 zeros, 0 poles\n      MISC: …\n     NOTES: …\n         T: …\n       FFT: …\nA non-empty FFTData looks like this:julia> julia> FFT\nFFTData with 188 ffts\n      NAME: \"TA.V04C..BHZ\"                     \n        ID: \"2006-02-01\"                       \n       LOC: 0.0 N, 0.0 E, 0.0 m\n        FS: 40.0\n      GAIN: 1.0\n   FREQMIN: 0.1\n   FREQMAX: 0.2\n    CC_LEN: 1800                               \n   CC_STEP: 450                                \n  WHITENED: false                              \n TIME_NORM: false                              \n      RESP: c = 0.0, 0 zeros, 0 poles\n      MISC: 0 entries                          \n     NOTES: 2 entries                          \n         T: 2006-02-01T00:07:30.000            …\n       FFT: 36001×188 Array{Complex{Float32}…     To access the Fourier transformed stored in FFT just do FFT.fftjulia> FFT.fft\n36001×188 Array{Complex{Float32},2}:\n 1.55183e-11+0.0im          1.59162e-12+0.0im          …  -9.66338e-12+0.0im        \n -7.41018e-7-7.1289e-9im    -4.28328e-7-6.37775e-9im        4.01078e-7-2.87796e-9im\n  -7.4106e-7-1.42547e-8im   -4.28291e-7-1.27635e-8im         4.0114e-7-5.75133e-9im\n -7.41128e-7-2.13963e-8im   -4.28259e-7-1.91511e-8im        4.01203e-7-8.65502e-9im\n -7.41227e-7-2.85335e-8im   -4.28155e-7-2.55307e-8im        4.01347e-7-1.15261e-8im\n -7.41355e-7-3.5675e-8im    -4.28056e-7-3.1921e-8im    …    4.01489e-7-1.44173e-8im\n -7.41517e-7-4.28364e-8im   -4.27914e-7-3.83248e-8im        4.01703e-7-1.7313e-8im  \n -7.41701e-7-5.00165e-8im   -4.27761e-7-4.47127e-8im        4.01912e-7-2.02371e-8im\n -7.41901e-7-5.72139e-8im   -4.27575e-7-5.11443e-8im        4.02202e-7-2.31354e-8im\n -7.42144e-7-6.44257e-8im   -4.27389e-7-5.75509e-8im        4.02492e-7-2.60687e-8im\n -7.42416e-7-7.16655e-8im   -4.27156e-7-6.39878e-8im   …    4.02823e-7-2.89802e-8im\n -7.42728e-7-7.89473e-8im   -4.26899e-7-7.04227e-8im        4.03202e-7-3.19224e-8im\n -7.43044e-7-8.62413e-8im   -4.26631e-7-7.68855e-8im        4.03602e-7-3.48864e-8im\n            ⋮                                          ⋱                            \n -3.96623e-9-2.01634e-11im  -7.62991e-9-1.17542e-11im      8.40682e-10-2.82507e-12im\n -3.96938e-9+1.15232e-11im  -7.63317e-9+1.61364e-11im      8.21519e-10-1.31251e-11im\n -3.97244e-9-7.90479e-13im  -7.63068e-9+1.09033e-11im  …   8.11657e-10-5.99587e-12im\n -3.96541e-9+1.46034e-11im  -7.61813e-9-6.24478e-12im      8.26425e-10-6.05316e-12im\n -3.97574e-9+2.11076e-12im  -7.64141e-9+7.56728e-12im      8.29292e-10+1.46063e-11im\n -3.95737e-9-1.91491e-12im  -7.62251e-9+1.36544e-11im       8.3114e-10+1.63336e-12im\n -3.97089e-9+2.54508e-12im  -7.63469e-9-1.97353e-12im      8.57792e-10+8.34888e-14im\n -3.97389e-9+6.95222e-13im  -7.62701e-9-1.36602e-11im  …    8.2856e-10+1.14622e-11im\n  -3.9615e-9+3.86602e-12im  -7.62552e-9+1.00606e-11im      8.41119e-10+1.14961e-11im\n  -3.9738e-9+2.04636e-12im  -7.64419e-9+1.50613e-12im       8.5253e-10-3.44774e-12im\n -3.96815e-9-8.12039e-12im  -7.63589e-9+2.64089e-12im       8.4199e-10-4.30272e-12im\n  -3.9762e-9+5.09676e-12im   -7.6402e-9-2.81547e-12im      8.31271e-10+3.22448e-12im\n -3.97375e-9+0.0im          -7.62793e-9+0.0im          …   8.40828e-10+0.0im      note: Note\nBy convention in Julia, data are stored column-wise. Here, FFT contains 188 individual Fourier transforms.FFTData"
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
    "text": "CorrData - Objects for ambient noise cross-correlations.Noise uses a custom structure CorrData for holding cross-correlations. To create an empty CorrData object, use the CorrData() function. CorrData are created with the compute_cc function.julia> using Noise\njulia> CorrData()\nCorrData with 0 Corrs\n      NAME: …\n        ID: …\n       LOC: 0.0 N, 0.0 E, 0.0 m\n      COMP: …\n   ROTATED: …\n CORR_TYPE: …\n        FS: 0.0\n      GAIN: 1.0\n   FREQMIN: 0.0\n   FREQMAX: 0.0\n    CC_LEN: …\n   CC_STEP: …\n  WHITENED: …\n TIME_NORM: …\n      RESP: c = 0.0, 0 zeros, 0 poles\n      MISC: …\n     NOTES: …\n    MAXLAG: 0.0\n         T: …\n      CORR: …A non-empty CorrData looks like this:julia> C\nCorrData with 188 Corrs\n      NAME: \"TA.V04C..BHZ.TA.V05C..BHZ\"        \n        ID: \"2006-02-01\"                       \n       LOC: 0.0 N, 0.0 E, 0.0 m\n      COMP: \"ZZ\"                               \n   ROTATED: false                              \n CORR_TYPE: \"coherence\"                        \n        FS: 40.0\n      GAIN: 1.0\n   FREQMIN: 0.1\n   FREQMAX: 0.2\n    CC_LEN: 1800                               \n   CC_STEP: 450                                \n  WHITENED: false                              \n TIME_NORM: false                              \n      RESP: c = 0.0, 0 zeros, 0 poles\n      MISC: 0 entries                          \n     NOTES: 2 entries                          \n    MAXLAG: 80.0\n         T: 2006-02-01T00:07:30.000            …\n      CORR: 6401×188 Array{Float32,2}     To access the correlations stored in C just do C.corrjulia> C.corr\n6401×188 Array{Float32,2}:\n 3.86871e-5   -0.00024922    0.000610002   0.000504283  …  -0.000187752   0.000239163   9.5009e-5  \n 3.8662e-5    -0.000238532   0.000607226   0.000493771     -0.000165946   0.000240476   7.08343e-5\n 3.81959e-5   -0.000227481   0.000602314   0.000481894     -0.000144085   0.000240927   4.63148e-5\n 3.76619e-5   -0.000216617   0.000596438   0.00046825      -0.000122814   0.000241393   2.16217e-5\n 3.6887e-5    -0.000205915   0.000588669   0.00045305      -0.000101321   0.000241209  -3.33542e-6\n 3.58405e-5   -0.000195269   0.000579826   0.000437172  …  -8.06184e-5    0.000240315  -2.85994e-5\n 3.4513e-5    -0.000184763   0.000569293   0.000419673     -5.94726e-5    0.00023961   -5.38829e-5\n 3.32238e-5   -0.000174677   0.000557559   0.000400775     -3.93204e-5    0.00023862   -7.91602e-5\n 3.15859e-5   -0.000164477   0.000543592   0.000381227     -1.84086e-5    0.000236609  -0.000104519\n 3.02724e-5   -0.00015467    0.000529009   0.000360099      1.71954e-6    0.000234192  -0.000129835\n 2.85284e-5   -0.00014533    0.000512517   0.000337998  …   2.13676e-5    0.000230932  -0.000155009\n 2.6588e-5    -0.000135612   0.000495403   0.000314866      4.03696e-5    0.000228341  -0.000180198\n 2.4749e-5    -0.000126298   0.000476726   0.000290888      5.9316e-5     0.000224295  -0.000205199\n ⋮                                                      ⋱   ⋮                                      \n 8.30166e-5    0.000492762   1.49007e-5   -0.000359526      0.000157452  -0.000232154  -0.000343711\n 8.6317e-5     0.00049216    5.2199e-6    -0.000358585      0.000152415  -0.000246062  -0.000354792\n 8.98544e-5    0.000490826  -4.95815e-6   -0.000357571  …   0.000146579  -0.000261209  -0.000365518\n 9.39397e-5    0.000489367  -1.49089e-5   -0.00035631       0.000140762  -0.000275435  -0.000375918\n 9.81829e-5    0.000486619  -2.55125e-5   -0.000354263      0.000135029  -0.000290028  -0.000385996\n 0.000102882   0.000483545  -3.53862e-5   -0.000352179      0.000128825  -0.000304768  -0.000395605\n 0.000107735   0.000480358  -4.58485e-5   -0.000349237      0.000122449  -0.000318608  -0.000404866\n 0.000112649   0.000476104  -5.58984e-5   -0.000346363  …   0.000116377  -0.000331822  -0.000413702\n 0.000118011   0.000471311  -6.60669e-5   -0.000342526      0.000110101  -0.000346527  -0.000422143\n 0.000123389   0.000465821  -7.64809e-5   -0.000338501      0.000103648  -0.000360307  -0.000430103\n 0.000128818   0.000459877  -8.66791e-5   -0.000334324      9.68363e-5   -0.000374563  -0.000437496\n 0.000134595   0.00045348   -9.59272e-5   -0.000329577      9.07089e-5   -0.000386802  -0.000444421\n 0.00014004    0.000446596  -0.000105674  -0.00032479   …   8.39108e-5   -0.000400697  -0.000450801note: Note\nBy convention in Julia, data are stored column-wise. Here, C contains 188 individual correlations. CorrData"
},

{
    "location": "preprocessing/#",
    "page": "Pre-Processing",
    "title": "Pre-Processing",
    "category": "page",
    "text": ""
},

{
    "location": "preprocessing/#Pre-Processing-methods-for-cleaning-raw-noise-data.-1",
    "page": "Pre-Processing",
    "title": "Pre-Processing - methods for cleaning raw noise data.",
    "category": "section",
    "text": "The pre-processing functions get raw data ready for correlation. Many of these rely on SeisIO\'s processing. A typical workflow includes the following steps:Loading data\nFilling gaps\nDownsampling\nSlicing day-long traces into smaller windows\nDemeaning, detrending and taperingHere is an example workflow:"
},

{
    "location": "preprocessing/#Load-or-download-data-1",
    "page": "Pre-Processing",
    "title": "Load or download data",
    "category": "section",
    "text": "Data can be downloaded using SeisIO\'s get_data function.julia> using SeisIO, Noise\njulia> S = get_data(\"IRIS\",\"TA.V04C..BHZ\",s=\"2006-02-01\",t=\"2006-02-02\")\nSeisData with 1 channels (1 shown)\n    ID: TA.V04C..BHZ                       \n  NAME: TA.V04C..BHZ                       \n   LOC: 0.0 N, 0.0 E, 0.0 m                \n    FS: 40.0                               \n  GAIN: 1.0                                \n  RESP: c = 0.0, 0 zeros, 0 poles          \n UNITS:                                    \n   SRC:                                    \n  MISC: 0 entries                          \n NOTES: 0 entries                          \n     T: 2006-02-01T00:00:00.020 (0 gaps)   \n     X: +1.044e-06                         \n        +1.060e-06                         \n            ...                            \n        +6.282e-07                         \n        (nx = 3456000)                     \n     C: 0 open, 0 total\njulia> writesac(S)Data can be read from disk using the read_data function .julia> S = read_data(\"sac\",\"2006.032.00.00.00.002.TA.V04C..BHZ.R.SAC\")\nSeisData with 1 channels (1 shown)\n    ID: TA.V04C..BHZ                       \n  NAME:                                    \n   LOC: 0.0 N, 0.0 E, 0.0 m                \n    FS: 40.0                               \n  GAIN: 1.0                                \n  RESP: c = 0.0, 0 zeros, 0 poles          \n UNITS:                                    \n   SRC: 2006.032.00.00.00.002.TA.V04C..BH…\n  MISC: 0 entries                          \n NOTES: 1 entries                          \n     T: 2006-02-01T00:00:00.002 (0 gaps)   \n     X: +1.044e-06                         \n        +1.060e-06                         \n            ...                            \n        +6.282e-07                         \n        (nx = 3456000)                     \n     C: 0 open, 0 total"
},

{
    "location": "preprocessing/#Merge-and-ungap-data-1",
    "page": "Pre-Processing",
    "title": "Merge and ungap data",
    "category": "section",
    "text": "Use the merge! method to merge S to have one channel. The ungap method replaces gaps in S[1].x with the mean of S[1].x.julia> merge!(S)\njulia> ungap!(S)"
},

{
    "location": "preprocessing/#Downsample-Data-1",
    "page": "Pre-Processing",
    "title": "Downsample Data",
    "category": "section",
    "text": "Downsample S to 20 Hz using process_raw!. process_raw!:Removes mean from each channel in S.\nDetrends each channel in S.\nDownsamples data to sampling rate fs\nPhase-shifts data to begin at 00:00:00.0\nThis is important for data that does not begin exactly on the sampling rate, e.g. the starttime is 2006-02-01T00:00:00.020.  julia> process_raw!(S,20.)\njulia> S\nSeisData with 1 channels (1 shown)\n    ID: TA.V04C..BHZ                       \n  NAME: TA.V04C..BHZ                       \n   LOC: 0.0 N, 0.0 E, 0.0 m                \n    FS: 20.0                               \n  GAIN: 1.0                                \n  RESP: c = 0.0, 0 zeros, 0 poles          \n UNITS:                                    \n   SRC:                                    \n  MISC: 0 entries                          \n NOTES: 1 entries                          \n     T: 2006-02-01T00:00:00.000 (0 gaps)   \n     X: +9.289e-10                         \n        -9.248e-10                         \n            ...                            \n        -9.287e-10                         \n        (nx = 1728000)                     \n     C: 0 open, 0 totalNote that now S has sampling rate fs = 20.0, has half as many points (nx = 1728001) as before and the starttime has changed to 2006-02-01T00:00:00.000."
},

{
    "location": "preprocessing/#Creating-Sliding-Windows-1",
    "page": "Pre-Processing",
    "title": "Creating Sliding Windows",
    "category": "section",
    "text": "Seats et al., 2011 showed that short term stacks of cross-correlations converge more quickly when dividing raw data into short, overlapping windows. The slide function splits S into windows of length cc_len (in seconds), with step between windows cc_step (in seconds).julia> cc_step, cc_len = 450, 1800 # 30 minute window length with 75% overlap\njulia> A, starts, ends = slide(S[1], cc_len, cc_step)\njulia> A\n36000×189 Array{Float64,2}:\n    0.0          -103.095  -67.9792   35.5011   -109.735   …   290.23     18.2172    2827.49    -8198.7      \n    3.14758e-6   -105.382  -69.2759   29.9752   -105.452       483.704    80.7089   -2596.52    -9547.97     \n    0.000105115  -104.929  -71.0974   24.1595   -104.596       635.558     9.51002  -1728.5     -7885.31     \n    0.000551506  -104.545  -73.0085   16.088    -100.576      1193.93     30.3235   -4105.61    -3688.67     \n    0.00155024   -104.646  -75.5247    6.76857   -95.9664     1181.32    -13.5895    1410.62      -65.8756   \n    ⋮                                                      ⋱     ⋮                                           \n -104.322        -105.706  -54.9128  -65.9102     51.2037  …   733.128  -422.694      -76.919      -5.99003  \n -106.025        -106.238  -63.692   -44.0787     46.1554      155.186  -378.362      -63.3664     -1.61041  \n -108.454        -106.888  -70.0719  -26.3757     42.7611     -367.41   -310.284     -142.125       0.503033\n -113.258        -104.777  -79.5973  -20.8923     38.9552      180.431  -272.636     -137.636       0.167523\n -113.7          -103.16   -85.4747   -7.70531    36.072       472.996  -293.358      -74.2375     -0.0133164Now A is a 2d Array containing 189 sliding windows, each 1800 seconds long."
},

{
    "location": "preprocessing/#Detrending-and-Demeaning-1",
    "page": "Pre-Processing",
    "title": "Detrending and Demeaning",
    "category": "section",
    "text": "The demean and detrend functions are applied column wise.julia> A = reshape(collect(1.:10.),5,2)\n5×2 Array{Float64,2}:\n 1.0   6.0\n 2.0   7.0\n 3.0   8.0\n 4.0   9.0\n 5.0  10.0\njulia> demean(A)\n5×2 Array{Float64,2}:\n -2.0  -2.0\n -1.0  -1.0\n  0.0   0.0\n  1.0   1.0\n  2.0   2.0\njulia> detrend(A)\n5×2 Array{Float64,2}:\n6.66134e-16  3.55271e-15\n1.33227e-15  5.32907e-15\n2.22045e-15  7.10543e-15\n2.66454e-15  7.10543e-15\n3.55271e-15  1.06581e-14"
},

{
    "location": "preprocessing/#Noise.process_raw!",
    "page": "Pre-Processing",
    "title": "Noise.process_raw!",
    "category": "function",
    "text": "process_raw!(S,fs)\n\nPre-process raw seismic data.\n\nRemoves mean from each channel in S.\nDetrends each channel in S.\nDownsamples data to sampling rate fs\nPhase-shifts data to begin at 00:00:00.0\n\nArguments\n\nS::SeisData: SeisData structure.\nfs::Float64: Sampling rate to downsample S.\n\n\n\n\n\n"
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
    "text": "slide(C::SeisChannel, cc_len::Int, cc_step::Int)\n\nCut C into equal length sliding windows.\n\nArguments\n\nC::SeisChannel: SeisChannel.\ncc_len::Int: Cross-correlation window length [s].\ncc_step::Int: Step between cross-correlation windows [s].\n\nReturns\n\nA::Array: Array of sliding windows\nstarts::Array: Array of start times of each window, in Unix time. E.g to convert       Unix time to date time, use u2d(starts[1]) = 2018-08-12T00:00:00\nends::Array: Array of end times of each window, in Unix time.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Noise.detrend!",
    "page": "Pre-Processing",
    "title": "Noise.detrend!",
    "category": "function",
    "text": "detrend!(X::AbstractArray{<:Union{Float32,Float64},1})\n\nRemove linear trend from array X using least-squares regression.\n\n\n\n\n\ndetrend!(X::AbstractArray{<:Union{Float32,Float64},2})\n\nRemove linear trend from columns of X using least-squares regression.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Noise.demean!",
    "page": "Pre-Processing",
    "title": "Noise.demean!",
    "category": "function",
    "text": "demean!(A::AbstractArray{<:Union{Float32,Float64},1})\n\nRemove mean from array A.\n\n\n\n\n\ndemean!(A::AbstractArray{<:Union{Float32,Float64},2})\n\nRemove mean from columns of array A.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Pre-Processing-Functions-1",
    "page": "Pre-Processing",
    "title": "Pre-Processing Functions",
    "category": "section",
    "text": "process_raw!\ndownsample\nslide\ndetrend!\ndemean!"
},

{
    "location": "filter/#",
    "page": "Filtering",
    "title": "Filtering",
    "category": "page",
    "text": ""
},

{
    "location": "filter/#Noise.bandpass!",
    "page": "Filtering",
    "title": "Noise.bandpass!",
    "category": "function",
    "text": "bandpass!(A,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandpass Filter.\n\nFilter data A from freqmin to freqmax using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreqmin::Float64: Pass band low corner frequency.\nfreqmax::Float64: Pass band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nbandpass!(C,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandpass Filter.\n\nFilter data in C from freqmin to freqmax using corners corners.\n\nArguments\n\nC::SeisChannel: SeisChannel to filter.\nfreqmin::Float64: Pass band low corner frequency.\nfreqmax::Float64: Pass band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nbandpass!(S,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandpass Filter.\n\nFilter channels in S from freqmin to freqmax using corners corners.\n\nArguments\n\nS::SeisData: SeisData to filter.\nfreqmin::Float64: Pass band low corner frequency.\nfreqmax::Float64: Pass band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "filter/#Noise.bandstop!",
    "page": "Filtering",
    "title": "Noise.bandstop!",
    "category": "function",
    "text": "bandstop!(A,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandstop Filter.\n\nFilter data A removing data between frequencies freqmin to freqmax using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreqmin::Float64: Stop band low corner frequency.\nfreqmax::Float64: Stop band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nbandstop!(C,freqmin,freqmax,fs,corners=4,zerophase=false)\n\nButterworth-Bandstop Filter.\n\nFilter data in C removing data between frequencies freqmin to freqmax using corners corners.\n\nArguments\n\nC::SeisChannel: SeisChannel to filter.\nfreqmin::Float64: Stop band low corner frequency.\nfreqmax::Float64: Stop band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "filter/#Noise.lowpass!",
    "page": "Filtering",
    "title": "Noise.lowpass!",
    "category": "function",
    "text": "lowpass(A,freq,fs,corners=4,zerophase=false)\n\nButterworth-Lowpass Filter.\n\nFilter data A over certain frequency freq using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nlowpass(C,freq,fs,corners=4,zerophase=false)\n\nButterworth-Lowpass Filter.\n\nFilter data in C over certain frequency freq using corners corners.\n\nArguments\n\nC::SeisChannel: SeisChannel to filter.\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nlowpass(S,freq,fs,corners=4,zerophase=false)\n\nButterworth-Lowpass Filter.\n\nFilter channels in S over certain frequency freq using corners corners.\n\nArguments\n\nS::SeisData: SeisData to filter.\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "filter/#Noise.highpass!",
    "page": "Filtering",
    "title": "Noise.highpass!",
    "category": "function",
    "text": "highpass(A,freq,fs,corners=4,zerophase=false)\n\nButterworth-Highpass Filter.\n\nFilter data A removing data below certain frequency freq using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nhighpass(C,freq,fs,corners=4,zerophase=false)\n\nButterworth-Highpass Filter.\n\nFilter data in C removing data below certain frequency freq using corners corners.\n\nArguments\n\nC::SeisChannel: SeisChannel to filter.\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\nhighpass(S,freq,fs,corners=4,zerophase=false)\n\nButterworth-Highpass Filter.\n\nFilter channels in S removing data below certain frequency freq using corners corners.\n\nArguments\n\nS::SeisData: SeisData to filter.\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "filter/#Filters-methods-for-filtering-SeisData-and-SeisChannel-objects.-1",
    "page": "Filtering",
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
    "location": "fft/#Noise.compute_fft",
    "page": "Computing FFTs",
    "title": "Noise.compute_fft",
    "category": "function",
    "text": "compute_fft(S::SeisData,freqmin::Float64,freqmax::Float64, fs::Float64\n                     cc_step::Int, cc_len::Int;\n                     time_norm::Union{Bool,String}=false,\n                     to_whiten::Bool=false)\n\nComputes windowed fft of ambient noise data.\n\nCross-correlates data using either cross-correlation, deconvolution or cross-coherence. Saves cross-correlations in JLD2 data set.\n\nArguments\n\nS::SeisChannel: SeisData structure.\nfreqmin::Float64: minimun frequency for whitening.\nfreqmax::Float64: maximum frequency for whitening.\nfs::Float64: Sampling rate to downsample S.\ncc_step::Int: time, in seconds, between successive cross-correlation windows.\ncc_len::Int: length of noise data window, in seconds, to cross-correlate.\ntime_norm::Union{Bool,String}: time domain normalization to perform.\nto_whiten::Bool: Apply whitening in frequency domain.\n\n\n\n\n\n"
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
    "text": "whiten\nprocess_fft\ncompute_fft\nsave_fft"
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
    "text": "clean_up!\ncorrelate\ncompute_cc\nsave_corr\nload_corr\nstack!"
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
    "text": "mwcs(ref, cur, fmin, fmax, fs, tmin, windowlength, windowstep,        smoothinghalfwin)\n\nChange in velocity measurement using the Moving Window Cross-Spectrum technique.\n\nThe current correlation cur is compared to the reference correlation ref. Both time series are sliced in several overlapping windows. Each slice is mean-adjusted and cosine-tapered (85% taper) before being Fourier- transformed to the frequency domain. F_cur(ν) and F_ref(ν) are the first halves of the Hermitian symmetric Fourier-transformed segments. The cross-spectrum X(ν) is defined as X(ν) = F_ref(ν) F_cur^*(ν) in which ^* denotes the complex conjugation. X(ν) is then smoothed by convolution with a Hanning window. The similarity of the two time-series is assessed using the cross-coherency between energy densities in the frequency domain: C(ν) = fracoverlineX(ν))sqrtoverlineF_ref(ν)^2 overlineF_cur(ν)^2 in which the over-line here represents the smoothing of the energy spectra for F_ref and F_cur and of the spectrum of X. The mean coherence for the segment is defined as the mean of C(ν) in the frequency range of interest. The time-delay between the two cross correlations is found in the unwrapped phase, ϕ(ν), of the cross spectrum and is linearly proportional to frequency: ϕ_j = m ν_j m = 2 π δ t The time shift for each window between two signals is the slope m of a weighted linear regression of the samples within the frequency band of interest. The weights are those introduced by [Clarke2011], which incorporate both the cross-spectral amplitude and cross-coherence, unlike [Poupinet1984]. The errors are estimated using the weights (thus the coherence) and the squared misfit to the modelled slope: e_m = sqrtsum_j(fracw_j ν_jsum_iw_i ν_i^2)^2σ_ϕ^2 where w are weights, ν are cross-coherences and σ_ϕ^2 is the squared misfit of the data to the modelled slope and is calculated as σ_ϕ^2 = fracsum_j(ϕ_j - m ν_j)^2N-1 The output of this process is a table containing, for each moving window: the central time lag, the measured delay, its error and the mean coherence of the segment.\n\nArguments\n\nref::AbstractArray: Reference correlation.\ncur::AbstractArray: Current correlation.\nfmin:Float64: minimum frequency in the correlation [Hz]\nfmax:Float64: maximum frequency in the correlation [Hz]\nfs:Float64: Sampling rate of ref and cur [Hz]\ntmin:Float64: The leftmost time lag [s]\nwindow_len:Float64: The moving window length [s]\nwindow_step:Float64: The step to jump for the moving window [s]\nsmoothing_half_win:Int: Defines the half length of the smoothing hanning                           window.\n\nReturns\n\ntime_axis:Array{Float64,1}: Central time of each moving window [s]\ndt:Array{Float64,1}: dt for each moving window\nerr:Array{Float64,1}: Errors for each moving window\nmcoh:Array{Float64,1}: Mean coherence for each moving window\n\nThis code is a Julia translation of the Python code from MSNoise.\n\n\n\n\n\n"
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

]}
