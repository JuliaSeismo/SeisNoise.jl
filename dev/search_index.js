var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#SeisNoise.jl-1",
    "page": "Home",
    "title": "SeisNoise.jl",
    "category": "section",
    "text": "Ambient Noise Cross-Correlation in Julia.SeisNoise.jl provides routines for quickly and efficiently implementing seismic interferometry.note: Note\nMuch of the I/O and preprocessing in SeisNoise.jl uses the SeisIO.jl package. Please read through the SeisIO Documentation to get familiar with seismic data processing in Julia."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "You can install the latest version of SeisNoise using the Julia package manager (Press ] to enter pkg). From the Julia command prompt:julia>]\n(@v1.4) pkg> add SeisNoiseOr, equivalently, via the Pkg API:julia> import Pkg; Pkg.add(\"SeisNoise\")We recommend using the latest version of SeisNoise by updating with the Julia package manager:(@v1.4) pkg> update SeisNoise"
},

{
    "location": "#Package-Features-1",
    "page": "Home",
    "title": "Package Features",
    "category": "section",
    "text": "Built upon SeisIO for easy and fast I/O.\nCustom structures for storing Raw Data, Fourier Transforms of data, and cross-correlations\nCPU/GPU compatible functions for cross-correlation.\nMethods for dv/v measurements.\nComing soon: Dispersion analysis."
},

{
    "location": "#Getting-Started-1",
    "page": "Home",
    "title": "Getting Started",
    "category": "section",
    "text": "SeisNoise.jl was designed to be as easy to use in the REPL as on an HPC cluster. This documentation provides a reference to all the underlying function for cross-correlation. We encourage you to develop your own workflow using SeisNoise\'s core functionality.(Image: plot1)"
},

{
    "location": "types/#",
    "page": "Workflow",
    "title": "Workflow",
    "category": "page",
    "text": ""
},

{
    "location": "types/#SeisNoise-Workflow-1",
    "page": "Workflow",
    "title": "SeisNoise Workflow",
    "category": "section",
    "text": "SeisNoise is designed around array-based cross-correlation. SeisNoise uses three custom data structures: RawData for storing windows ambient noise, FFTData for storing Fourier transforms of ambient noise, and CorrData for storing cross-correlations. SeisNoise types extend SeisIO\'s SeisChannel type for 2D-array ambient noise processing.(Image: dataflow)note: Note\nBy convention in Julia, array data is stored column-wise. Noise windows in   RawData objects, ffts in FFTData objects, and cross-corrleations in   CorrData objects are all stored column-wise. To access data one or more   time windows, use column-based indexing, e.g. R[:,1], F[:,2], or C[:,3:5]."
},

{
    "location": "types/#Timing-1",
    "page": "Workflow",
    "title": "Timing",
    "category": "section",
    "text": "SeisNoise stores the start time of each window stored in RawData, FFTData and CorrData objects in the .t field as a 1D array. Timestamps are stored according to Unix Time in double-precision floating-point format (Float64). SeisNoise uses the Dates module for all timing. For instance, January 1st, 2020 in Unix time is given byusing Dates\ndatetime2unix(DateTime(2020,1,1))\n1.5778368e9where datetime2unix is a function in the Dates module for converting between Unix Time and DateTime structures. Dates also provides a unix2datetime for converting from Unix Time to DateTime structures. Let\'s say we have an array of timestamps, t, with hourly samples between 01/01/2020 and 01/02/2020:t = 1.5778368e9 .+ collect(0:3600:3600*24)\n25-element Array{Float64,1}:\n 1.5778368e9\n 1.5778404e9\n 1.577844e9\n 1.5778476e9\n 1.5778512e9\n 1.5778548e9\n ⋮\n 1.5779088e9\n 1.5779124e9\n 1.577916e9\n 1.5779196e9\n 1.5779232e9to convert to timestamps, we would use the unix2datetime function along with the Julia broadcast operator, ., as shown below:timestamps = unix2datetime.(t)\n25-element Array{DateTime,1}:\n 2020-01-01T00:00:00\n 2020-01-01T01:00:00\n 2020-01-01T02:00:00\n 2020-01-01T03:00:00\n 2020-01-01T04:00:00\n 2020-01-01T05:00:00\n ⋮\n 2020-01-01T20:00:00\n 2020-01-01T21:00:00\n 2020-01-01T22:00:00\n 2020-01-01T23:00:00\n 2020-01-02T00:00:00"
},

{
    "location": "types/#RawData-Objects-for-windowed-raw-data-1",
    "page": "Workflow",
    "title": "RawData - Objects for windowed raw data",
    "category": "section",
    "text": "RawData objects store windows of ambient noise. RawData have very similar structure to SeisChannel objects, except with added fields for:cc_len: the cross-correlation window length (in seconds).\ncc_step: the step length between successive correlation windows (in seconds).\nfreqmin: minimum frequency RawData has been filtered (in Hz).\nfreqmax: maximum frequency RawData has been filtered (in Hz).\ntime_norm: Type of amplitude normalization applied (onebit, clipping, etc..).\nt: Starttime of each noise window (number of seconds since the unix epoch 1970-01-01T00:00:00 as a Float64).\nx: Raw data stored in columns (2d-array).To create an empty RawData object, use the RawData() function. Here is an empty RawData object in the REPL:julia> using SeisNoise\njulia> RawData()\nRawData with 0 windows\n      NAME: …\n        ID: …\n       LOC: 0.0 N, 0.0 E, 0.0 m\n        FS: 0.0\n      GAIN: 1.0\n   FREQMIN: 0.0\n   FREQMAX: 0.0\n    CC_LEN: …\n   CC_STEP: …\n TIME_NORM: …\n      RESP: a0 1.0, f0 1.0, 0z, 0p\n      MISC: …\n     NOTES: …\n         T: …\n         X: …A more natural way to create a RawDataobject is to convert an existing SeisChannel or SeisData object into a RawData object. This allows one to do preprocessing that is more natural for single component data (e.g. instrument response removal) before doing array-based processing. Here is an example of converting a SeisData object to a RawData object.note: Note\ncc_len and cc_step must be included to create a RawData object.julia> using SeisIO, SeisNoise\n\njulia> S = get_data(\"FDSN\", \"TA.M22K..BHZ\", src=\"IRIS\",s=\"2019-01-01\",t=600)\nSeisData with 1 channels (1 shown)\n    ID: TA.M22K..BHZ                       \n  NAME: Willow, AK, USA                    \n   LOC: 61.7531 N, -150.12 E, 57.0 m       \n    FS: 40.0                               \n  GAIN: 5.03883e8                          \n  RESP: a0 8.32666e17, f0 0.2, 6z, 11p     \n UNITS: m/s                                \n   SRC: http://service.iris.edu/fdsnws/da…\n  MISC: 4 entries                          \n NOTES: 0 entries                          \n     T: 2019-01-01T00:00:00.000 (0 gaps)   \n     X: +8.540e+02                         \n        +6.230e+02                         \n            ...                            \n        +1.081e+03                         \n        (nx = 24001)                       \n     C: 0 open, 0 total\n\njulia> cc_len, cc_step = 100.,50. # specify window length and step (in seconds)\n(100, 50)\n\njulia> R = RawData(S,cc_len,cc_step)\nRawData with 11 windows\n      NAME: \"Willow, AK, USA\"                  \n        ID: \"TA.M22K..BHZ\"                     \n       LOC: 61.7531 N, -150.12 E, 57.0 m\n        FS: 40.0\n      GAIN: 5.03883e8\n   FREQMIN: 0.01\n   FREQMAX: 20.0\n    CC_LEN: 100.                                \n   CC_STEP: 50.                                 \n TIME_NORM: false                              \n      RESP: a0 8.32666e17, f0 0.2, 6z, 11p\n      MISC: 4 entries                          \n     NOTES: 2 entries                          \n         T: 2019-01-01T00:00:00.000            …\n         X: 4000×11 Array{Float32,2}     \nThis created a RawData object, R, 11 windows, each 100 seconds long and sampled at 40 Hz (4000 points). The freqmin is automatically set to 0.01 Hz because this is the 1 / cc_len, whereas freqmax is set to the Nyquist frequency R.fs/2. To access the noise data stored in R just do R.xjulia> R.x\n4000×11 Array{Float32,2}:\n  854.0  1645.0   790.0  976.0   959.0  …   921.0  1148.0   718.0\n  623.0   292.0   762.0  790.0  1065.0      828.0  1242.0   596.0\n  530.0  1398.0   408.0  649.0  1033.0      -32.0  1263.0   678.0\n  684.0   970.0   546.0  793.0   984.0     1399.0  1324.0   676.0\n  810.0   273.0  1066.0  943.0   992.0      561.0  1340.0   655.0\n  834.0   554.0  1385.0  866.0   955.0  …    62.0  1291.0   729.0\n    ⋮                                   ⋱                     ⋮  \n  647.0   810.0   801.0  380.0  1683.0  …   830.0   815.0  1223.0\n  639.0   700.0   880.0  285.0  1284.0      709.0   798.0  1149.0\n  860.0   789.0   955.0  133.0  1427.0      844.0   913.0   998.0\n 1052.0   875.0   927.0  115.0  1623.0      857.0   944.0   970.0\n  872.0   925.0   883.0    2.0  1141.0      810.0   883.0  1061.0See the Pre-Processing and Filtering pages for functions that work on RawData objects.  "
},

{
    "location": "types/#FFTData-Objects-for-Fourier-transforms-(FFTs)-1",
    "page": "Workflow",
    "title": "FFTData - Objects for Fourier transforms (FFTs)",
    "category": "section",
    "text": "Cross-correlation in SeisNoise is done in the frequency domain:C_AB(ω) = u^*_A(ω) u_B(ω)where C_AB(ω) is the element-wise multiplication between u^*_A(ω), the complex conjugate of the Fourier transform of ambient noise time series A, and u_B(ω), the Fourier transform of ambient noise time series B. This exploits the O(nlogn) complexity of the Fast-Fourier Transforms. Since ambient noise data is real valued, SeisNoise uses real Fast-Fourier Transforms (rfft), which offer a 2-3x speed over a general fft.The FFTData structure in SeisNoise stores Fourier transforms (u(ω)) of ambient noise. FFTData allow users to apply smoothing operations, such as whitening, coherence, or deconvolution, in-place before cross-corrleating. To create an empty FFTData object, use the FFTData() function.julia> using SeisNoise\njulia> FFTData()\nFFTData with 0 ffts\n      NAME: …\n        ID: …\n       LOC: 0.0 N, 0.0 E, 0.0 m\n        FS: 0.0\n      GAIN: 1.0\n   FREQMIN: 0.0\n   FREQMAX: 0.0\n    CC_LEN: …\n   CC_STEP: …\n  WHITENED: …\n TIME_NORM: …\n      RESP: c = 0.0, 0 zeros, 0 poles\n      MISC: …\n     NOTES: …\n         T: …\n       FFT: …\nThe only difference between an FFTData and a RawData is the swapping of the x data field to the fft data field.FFTData more naturally flow from RawData input to the  rfft function. julia> F = rfft(R)\nFFTData with 11 ffts\n      NAME: \"TA.M22K..BHZ\"                     \n        ID: \"2019-01-01\"                       \n       LOC: 61.7531 N, -150.12 E, 57.0 m\n        FS: 40.0\n      GAIN: 5.03883e8\n   FREQMIN: 0.01\n   FREQMAX: 20.0\n    CC_LEN: 100.                                \n   CC_STEP: 50.                                 \n  WHITENED: false                              \n TIME_NORM: false                              \n      RESP: a0 8.32666e17, f0 0.2, 6z, 11p\n      MISC: 4 entries                          \n     NOTES: 2 entries                          \n         T: 2019-01-01T00:00:00.000            …\n       FFT: 2001×11 Array{Complex{Float32},2}  To access the Fourier transform stored in F just do F.fftjulia> F.fft\n2001×11 Array{Complex{Float32},2}:\n 3.18286e6+0.0im      3.18685e6+0.0im       …  3.19561e6+0.0im    \n  -3205.34-7445.99im    12079.2+1573.37im        11842.9-4106.75im\n  -10263.9+5427.09im    -1056.6+3629.68im        5421.51+2751.05im\n  -6255.81-25112.9im    21741.7+9154.63im        37171.6+10855.5im\n  -6164.36-7888.87im    27625.1-2106.65im        10308.3-24606.1im\n   2470.13+8848.62im    13407.9+5940.08im   …   -10548.2+1305.7im\n          ⋮                                 ⋱           ⋮         \n   133.262-18.7686im    1204.43+18.5029im       -135.982+7.42407im\n   122.472-7.24219im    1179.84+1.39209im       -133.742-5.74707im\n   131.036+7.63867im    1199.29+0.913574im      -117.602+4.26135im\n   123.142+4.82568im    1202.13+1.8291im        -112.404+11.0547im\n     126.0+0.0im         1179.0+0.0im       …     -148.0+0.0im    "
},

{
    "location": "types/#CorrData-Objects-for-ambient-noise-cross-correlations-1",
    "page": "Workflow",
    "title": "CorrData - Objects for ambient noise cross-correlations",
    "category": "section",
    "text": "As stated in the previous section, cross-correlation in the frequency domain is an element-wise product of Fourier-transforms of ambient noise:C_AB(ω) = u^*_A(ω) u_B(ω)To transform a frequency domain cross-correlation into the time domain, we take the inverse real Fourier transform of C_AB(ω):C(τ)_AB = mathfrakF^-1 left(C_AB(ω)right)where τ is the lag time. Computing a cross-correlation thus necessitates two Fourier transforms, one element-wise multiplication (of complex numbers), and an inverse Fourier transform.The CorrData structure in SeisNoise stores time-domain cross-correlations. To create an empty CorrData object, use the CorrData() function.julia> CorrData()\nCorrData with 0 Corrs\n      NAME: …\n        ID: …\n       LOC: 0.0 N, 0.0 E, 0.0 m\n      COMP: …\n   ROTATED: …\n CORR_TYPE: …\n        FS: 0.0\n      GAIN: 1.0\n   FREQMIN: 0.0\n   FREQMAX: 0.0\n    CC_LEN: …\n   CC_STEP: …\n  WHITENED: …\n TIME_NORM: …\n      RESP: a0 1.0, f0 1.0, 0z, 0p\n      MISC: …\n     NOTES: …\n      DIST: 0.0\n       AZI: 0.0\n       BAZ: 0.0\n    MAXLAG: 0.0\n         T: …\n      CORR: …CorrData have a few additional parameters:comp: The component of the cross-correlations, e.g. EN, RT, ZZ, etc...\nrotated: True/false if correlation has been rotated.\ncorr_type: The type of correlation, e.g. \"cross-correlation\", \"coherence\", or \"deconvolution\".\ndist: Distance between station 1 and station 2 (in Km).\nazi: Azimuth from station 1 and station 2 (in degrees).\nbaz: Back azimuth between station 1 and station 2 (in degrees).\nmaxlag: Maximum lag time (in seconds) in cross-correlations.\ncorr: Time-domain cross-correlations stored in columns.CorrData more naturally flow from FFTData input to the compute_cc function. In this example, we are computing an \"autocorrleation\" by inputting the same an FFTData twice:julia> maxlag = 20.\n20.0\n\njulia> C = correlate(F,F,maxlag,corr_type=\"cross-correlation\")\nCorrData with 11 Corrs\n      NAME: \"TA.M22K..BHZ.TA.M22K..BHZ\"        \n        ID: \"2019-01-01\"                       \n       LOC: 61.7531 N, -150.12 E, 57.0 m\n      COMP: \"ZZ\"                               \n   ROTATED: false                              \n CORR_TYPE: \"cross-correlation\"                \n        FS: 40.0\n      GAIN: 5.03883e8\n   FREQMIN: 0.01\n   FREQMAX: 20.0\n    CC_LEN: 100.                                \n   CC_STEP: 50.                                 \n  WHITENED: false                              \n TIME_NORM: false                              \n      RESP: a0 8.32666e17, f0 0.2, 6z, 11p\n      MISC: 4 entries                          \n     NOTES: 2 entries                          \n      DIST: 0.0\n       AZI: 0.0\n       BAZ: 0.0\n    MAXLAG: 20.0\n         T: 2019-01-01T00:00:00.000            …\n      CORR: 1601×11 Array{Float32,2}   To access the autocorrelations stored in C just do C.corr:julia> C.corr\n1601×11 Array{Float32,2}:\n 2.56154e9  2.5478e9   2.51531e9  …  2.53951e9  2.55251e9  2.51616e9\n 2.54715e9  2.52733e9  2.51668e9     2.54039e9  2.55533e9  2.51614e9\n 2.53042e9  2.54303e9  2.52005e9     2.53713e9  2.55813e9  2.51694e9\n 2.57688e9  2.57981e9  2.52313e9     2.53812e9  2.55982e9  2.52244e9\n 2.617e9    2.57287e9  2.52175e9     2.53825e9  2.56081e9  2.52769e9\n 2.51533e9  2.50786e9  2.5218e9   …  2.53502e9  2.56186e9  2.52965e9\n 2.5034e9   2.52427e9  2.52484e9     2.53588e9  2.56324e9  2.53052e9\n 2.5998e9   2.578e9    2.52705e9     2.53676e9  2.56418e9  2.53246e9\n ⋮                                ⋱                        ⋮        \n 2.56922e9  2.56349e9  2.52909e9     2.53253e9  2.56514e9  2.53592e9\n 2.5998e9   2.578e9    2.52705e9  …  2.53676e9  2.56418e9  2.53246e9\n 2.5034e9   2.52427e9  2.52484e9     2.53588e9  2.56324e9  2.53052e9\n 2.51533e9  2.50786e9  2.5218e9      2.53502e9  2.56187e9  2.52965e9\n 2.617e9    2.57287e9  2.52175e9     2.53825e9  2.56081e9  2.52769e9\n 2.57688e9  2.57981e9  2.52313e9     2.53812e9  2.55982e9  2.52244e9\n 2.53042e9  2.54303e9  2.52005e9  …  2.53713e9  2.55813e9  2.51694e9"
},

{
    "location": "types/#SeisNoise-Functions-1",
    "page": "Workflow",
    "title": "SeisNoise Functions",
    "category": "section",
    "text": "SeisNoise\'s modular functions work on RawData, FFTData, and CorrData objects in-place through Julia\'s multiple-dispatch. Functions in SeisNoise that end in ! (e.g. onebit!(R)) by convention modify their arguments, while functions that do not end in ! (e.g. onebit(R)) allocate new arrays or objects. Have a look at the Pre-Processing, Computing FFTs and Cross-Correlating sections to learn how to apply functions to NoiseData objects."
},

{
    "location": "types/#SeisNoise.RawData",
    "page": "Workflow",
    "title": "SeisNoise.RawData",
    "category": "type",
    "text": "RawData\n\nA structure for raw ambient noise data.\n\nFields: RawData\n\nField Description\n:name Freeform channel names\n:id Channel ids. use NET.STA.LOC.CHAN format when possible.\n:loc Location (position) object\n:fs Sampling frequency in Hz.\n:gain Scalar gain; divide data by the gain to convert to units\n:freqmin Minimum frequency from instrument response.\n:freqmax Maximum frequency from instrument response.\n:cc_len Raw data window length in seconds.\n:cc_step Spacing between windows in seconds.\n:time_norm Apply one-bit whitening with \"one_bit\".\n:resp SeisIO InstrumentResponse object\n:misc Dictionary for non-critical information.\n:notes Timestamped notes; includes automatically-logged acquisition and\n processing information.\n:t Starttime of each window.\n:x Raw data stored in columns.\n\n\n\n\n\n"
},

{
    "location": "types/#SeisNoise.FFTData",
    "page": "Workflow",
    "title": "SeisNoise.FFTData",
    "category": "type",
    "text": "FFTData\n\nA structure for fourier transforms (FFT) of ambient noise data.\n\nFields: FFTData\n\nField Description\n:name Freeform channel names\n:id Channel ids. use NET.STA.LOC.CHAN format when possible.\n:loc Location (position) object\n:fs Sampling frequency in Hz.\n:gain Scalar gain; divide data by the gain to convert to units\n:freqmin Minimum frequency for whitening.\n:freqmax Maximum frequency for whitening.\n:cc_len Raw data window length in number of points.\n:cc_step Spacing between windows in number of points.\n:whitened Whitening applied.\n:time_norm Apply one-bit whitening with \"one_bit\".\n:resp SeisIO InstrumentResponse object\n:misc Dictionary for non-critical information.\n:notes Timestamped notes; includes automatically-logged acquisition and\n processing information.\n:t Starttime of each FFT.\n:fft FFTs stored in columns.\n\n\n\n\n\n"
},

{
    "location": "types/#SeisNoise.CorrData",
    "page": "Workflow",
    "title": "SeisNoise.CorrData",
    "category": "type",
    "text": "CorrData\n\nA structure for cross-correlations of ambient noise data.\n\nFields: CorrData\n\nField Description\n:name Freeform channel name in NET1.STA1.LOC1.CHAN1.NET2.STA2.LOC2.CHAN2 form\n:id Channel ids. use NET.STA.LOC.CHAN format when possible.\n:loc Location (position) object for station 1.\n:comp Component of cross-correlation (e.g ZZ, RT, etc..).\n:rotated Boolean for rotation.\n:fs Sampling frequency in Hz.\n:gain Scalar gain; divide data by the gain to convert to units.\n:freqmin Minimum frequency for whitening.\n:freqmax Maximum frequency for whitening.\n:cc_len Raw data window length in number of points.\n:cc_step Spacing between windows in number of points.\n:whitened Whitening applied.\n:time_norm Apply one-bit whitening with \"one_bit\".\n:resp SeisIO InstrumentResponse object\n:misc Dictionary for non-critical information.\n:notes Timestamped notes; includes automatically-logged acquisition and\n processing information.\n:dist Distance between station 1 and station 2 in Km.\n:azi Azimuth from station 1 and station 2 in degrees.\n:baz Back azimuth between station 1 and station 2 in degrees.\n:maxlag Maximum lag time in seconds to keep from correlations.\n:t Starttime of each correlation.\n:corr Correlations stored in columns.\n\n\n\n\n\n"
},

{
    "location": "types/#Type-Documentation-1",
    "page": "Workflow",
    "title": "Type Documentation",
    "category": "section",
    "text": "RawData\nFFTData\nCorrData"
},

{
    "location": "using_gpus/#",
    "page": "Using the GPU",
    "title": "Using the GPU",
    "category": "page",
    "text": ""
},

{
    "location": "using_gpus/#Using-GPUs-1",
    "page": "Using the GPU",
    "title": "Using GPUs",
    "category": "section",
    "text": "SeisNoise is written to run on graphical processing units (GPUs) for increased performance. Depending on your CPU and GPU combination, speedups of >10-20x are possible.tip: Running on GPUs\nIf you are having issues with running SeisNoise on a GPU or setting things up, please open an issue and we\'ll do our best to help out!"
},

{
    "location": "using_gpus/#When-to-use-a-GPU-1",
    "page": "Using the GPU",
    "title": "When to use a GPU",
    "category": "section",
    "text": "GPUs are very useful for number crunching. If you have a large number of stations to cross-correlate, using a GPU could greatly improve your cross-correlation processing time.GPU processing tends to be memory-limited. High-end GPUs (think $$) such as the Nvidia Tesla V100 only come with up to 32 GB of memory. On a GPU with 16 GB of memory, you can cross-correlate ~100 day-long 3-component stations sampled at 100 Hz in memory."
},

{
    "location": "using_gpus/#Getting-access-to-GPUs-1",
    "page": "Using the GPU",
    "title": "Getting access to GPUs",
    "category": "section",
    "text": "If you don\'t have a GPU there are a few resources you can try to acquire a GPU from.If you have access to any supercomputer clusters, check to see if they have any GPUs. See also this Stack Overflow post: Where can I get access to GPU cluster for educational purpose?Cloud computing providers such as Google Cloud and Amazon EC2 allow you to rent GPUs per hour. Sometimes they offer free trials or credits that can be used towards GPUs although they seem to be getting less common.See the Julia on Google Colab: Free GPU-Accelerated Shareable Notebooks post on the Julia Discourse.Code Ocean also has GPU support and allows you to spin up capsules with pretty decent Tesla K80 GPUs for free (for now) if you want to play around with them.  You\'ll want to use their \"Ubuntu Linux with GPU support (18.04.3)\" with the ability to compile CUDA code. Then you\'ll have to install Julia manually.NextJournal provides a Jupyter-like Julia inferface with access to Tesla K80 GPUs for free."
},

{
    "location": "using_gpus/#I-have-a-GPU.-Now-what?-1",
    "page": "Using the GPU",
    "title": "I have a GPU. Now what?",
    "category": "section",
    "text": "Make sure you have an Nvidia GPU that is CUDA compatible: https://developer.nvidia.com/cuda-gpus. Most recent GPUs should be but older GPUs and many laptop GPUs may not be.Then download and install the CUDA Toolkit: https://developer.nvidia.com/cuda-downloadsOnce the CUDA Toolkit is installed, you might have to build SeisNoise againjulia>]\n(v1.4) pkg> build SeisNoise"
},

{
    "location": "using_gpus/#Using-the-GPU-with-Julia-1",
    "page": "Using the GPU",
    "title": "Using the GPU with Julia",
    "category": "section",
    "text": "SeisNoise uses the Cuda.jl library to do all GPU-based operations. To get started with GPU-programming, take a look at the Cuda.jl introduction tutorial. Here is a great slide from Tim Besard on GPU computing in Julia.Allocating arrays on the GPU is simple using the cu function. This example creates a random maxtrix on the CPU, then sends it to the GPU using cu:using CUDA\nA = cu(rand(1000,10))\n1000×10 CuArray{Float32,2,Nothing}:\n 0.483598   0.0278827  0.564381  …  0.0179997  0.0894625  0.943952\n 0.532773   0.601974   0.965769     0.91777    0.746491   0.837131\n 0.0613844  0.638339   0.401549     0.385933   0.627327   0.682377\n 0.597944   0.349372   0.147358     0.607321   0.592268   0.381262\n 0.516212   0.505277   0.35825      0.0874267  0.946768   0.525376\n ⋮                               ⋱                        \n 0.696375   0.694293   0.646136     0.0734224  0.0999966  0.804346\n 0.266413   0.710535   0.762904     0.13076    0.0786623  0.914678\n 0.255818   0.0446712  0.943867     0.981916   0.835427   0.506067\n 0.611743   0.715961   0.181798     0.793934   0.568452   0.84297\nor similarly using the Julia pipe syntax,A = rand(1000,10) |> cu\n1000×10 CuArray{Float32,2,Nothing}:\n 0.511658   0.231007   0.819      …  0.110045   0.517574    0.195391\n 0.589017   0.69037    0.0288187     0.599881   0.0290229   0.136348\n 0.175985   0.217625   0.124626      0.0619627  0.862026    0.431731\n 0.0430858  0.569207   0.585091      0.605072   0.333407    0.505581\n 0.510674   0.969721   0.875608      0.984235   0.0604092   0.0494258\n ⋮                                ⋱                         \n 0.938167   0.673817   0.0642641     0.401075   0.420085    0.0400348\n 0.311241   0.0985884  0.309681      0.33953    0.793376    0.385521\n 0.738205   0.500696   0.294173      0.681627   0.00480331  0.751362\n 0.740553   0.566934   0.845259      0.503732   0.970519    0.218747\nNote that here A is a CuArray which is stored in memory on the GPU. Most standard Julia functions (matrix multiplication, FFT, mean, max, etc..) will work out of the box on CuArrays. For example, here is matrix multiplicationB = rand(10,1000) |> cu\nA * B\n1000×1000 CuArray{Float32,2,Nothing}:\n 3.02011  2.85182  3.17944  …  2.20066  1.89412  2.32329\n 3.12874  2.83144  3.34268     1.65363  1.97212  2.41745\n 3.22296  2.53652  2.84638     2.11542  2.31627  1.71027\n 3.21165  2.51493  2.9726      2.07334  2.55837  2.63839\n 3.2731   2.83646  3.75265     1.44419  1.68725  3.5038\n 2.18507  1.80968  2.13033  …  1.47744  1.97971  1.59087\n ⋮                          ⋱                    \n 1.95406  1.95351  2.62476  …  1.07703  1.49473  1.72594\n 2.26404  2.09106  2.86892     1.12002  1.21775  1.73231\n 2.56381  2.1512   2.70032     2.03604  1.98124  1.38784\n 3.07007  3.02482  3.36264     2.14246  2.4654   2.34386\n 3.46125  2.72799  3.75003     1.8953   2.44945  2.65501"
},

{
    "location": "using_gpus/#SeisNoise-on-the-GPU-1",
    "page": "Using the GPU",
    "title": "SeisNoise on the GPU",
    "category": "section",
    "text": "SeisNoise can process data and compute cross-correlations on the GPU with CUDA. CUDA.jl provides an the (CuArray) type for storing data on the GPU. Data in SeisNoise structures (R.x, F.fft, and C.corr fields, for RawData, FFTData, and CorrData, respectively) can move between an Array on the CPU to a CuArray on the GPU using the gpu and cpu functions, as shown below.   note: Note\nOnly Nvidia GPUs are suported at the moment. Hold in there for AMD/OpenCL/Metal support...# create raw data and send to GPU\nR = RawData(S1, cc_len, cc_step) |> gpu\nR.x\n72000×188 CuArrays.CuArray{Float32,2,Nothing}\n\n# send data back to the CPU\nR = R |> cpu\nR.x\n72000×188 Array{Float32,2}All basic processing remains the same on the GPU. Here is a complete cross-correlation routine on the GPU:# send data to GPU\nR1 = RawData(S1, cc_len, cc_step) |> gpu\nR2 = RawData(S2, cc_len, cc_step) |> gpu\nR = [R1,R2]\n\n# preprocess on the GPU\ndetrend!.(R)\ntaper!.(R)\nbandpass!.(R,freqmin,freqmax,zerophase=true)\n\n# Real FFT on GPU\nFFT = rfft.(R)\nwhiten!.(FFT,freqmin,freqmax)\n\n# compute correlation and send to cpu\nC = correlate(FFT[1],FFT[2],maxlag) |> cpu"
},

{
    "location": "using_gpus/#Routines-Implemented-on-the-GPU-1",
    "page": "Using the GPU",
    "title": "Routines Implemented on the GPU",
    "category": "section",
    "text": "Preprocessing:\ndetrend,demean, taper, onebit, smooth\nFiltering:\nbandpass, bandstop, lowpass, highpass\nFourier Domain:\nwhiten, rfft, irfft\nCross-correlation:\ncorrelate, cross-coherence, deconvolution\nPost-processing:\nstack, filters, etc.."
},

{
    "location": "using_gpus/#Credits-1",
    "page": "Using the GPU",
    "title": "Credits",
    "category": "section",
    "text": "This page builds heavily on the using GPUs page from Oceananigans."
},

{
    "location": "preprocessing/#",
    "page": "Pre-Processing",
    "title": "Pre-Processing",
    "category": "page",
    "text": ""
},

{
    "location": "preprocessing/#Pre-Processing-1",
    "page": "Pre-Processing",
    "title": "Pre-Processing",
    "category": "section",
    "text": "The pre-processing functions get raw data ready for correlation. Many of these rely on SeisIO\'s processing. A typical workflow includes the following steps:Loading data\nFilling gaps\nDownsampling\nSlicing day-long traces into smaller windows\nDemeaning, detrending and tapering\nTime-domain normalizationHere is an example workflow:"
},

{
    "location": "preprocessing/#Load-or-download-data-1",
    "page": "Pre-Processing",
    "title": "Load or download data",
    "category": "section",
    "text": "Data can be downloaded using SeisIO\'s get_data function.julia> using SeisIO, SeisNoise\njulia> S = get_data(\"IRIS\",\"TA.V04C..BHZ\",s=\"2006-02-01\",t=\"2006-02-02\")\nSeisData with 1 channels (1 shown)\n    ID: TA.V04C..BHZ                       \n  NAME: TA.V04C..BHZ                       \n   LOC: 0.0 N, 0.0 E, 0.0 m                \n    FS: 40.0                               \n  GAIN: 1.0                                \n  RESP: c = 0.0, 0 zeros, 0 poles          \n UNITS:                                    \n   SRC:                                    \n  MISC: 0 entries                          \n NOTES: 0 entries                          \n     T: 2006-02-01T00:00:00.020 (0 gaps)   \n     X: +1.044e-06                         \n        +1.060e-06                         \n            ...                            \n        +6.282e-07                         \n        (nx = 3456000)                     \n     C: 0 open, 0 total\njulia> writesac(S)Data can be read from disk using the read_data function .julia> S = read_data(\"sac\",\"2006.032.00.00.00.002.TA.V04C..BHZ.R.SAC\")\nSeisData with 1 channels (1 shown)\n    ID: TA.V04C..BHZ                       \n  NAME:                                    \n   LOC: 0.0 N, 0.0 E, 0.0 m                \n    FS: 40.0                               \n  GAIN: 1.0                                \n  RESP: c = 0.0, 0 zeros, 0 poles          \n UNITS:                                    \n   SRC: 2006.032.00.00.00.002.TA.V04C..BH…\n  MISC: 0 entries                          \n NOTES: 1 entries                          \n     T: 2006-02-01T00:00:00.002 (0 gaps)   \n     X: +1.044e-06                         \n        +1.060e-06                         \n            ...                            \n        +6.282e-07                         \n        (nx = 3456000)                     \n     C: 0 open, 0 total"
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
    "location": "preprocessing/#Detrending-and-Demeaning-1",
    "page": "Pre-Processing",
    "title": "Detrending and Demeaning",
    "category": "section",
    "text": "The demean and detrend functions are applied column wise.julia> A = reshape(collect(1.:10.),5,2)\n5×2 Array{Float64,2}:\n 1.0   6.0\n 2.0   7.0\n 3.0   8.0\n 4.0   9.0\n 5.0  10.0\njulia> demean(A)\n5×2 Array{Float64,2}:\n -2.0  -2.0\n -1.0  -1.0\n  0.0   0.0\n  1.0   1.0\n  2.0   2.0\njulia> detrend(A)\n5×2 Array{Float64,2}:\n6.66134e-16  3.55271e-15\n1.33227e-15  5.32907e-15\n2.22045e-15  7.10543e-15\n2.66454e-15  7.10543e-15\n3.55271e-15  1.06581e-14demean and detrend also work on RawData and CorrData structure. If R is a RawData, then the windowed data stored in R.x can be detrended in-place usingdetrend!(R)or allocated to a new RawData usingRd = detrend(R)"
},

{
    "location": "preprocessing/#Amplitude-Normaliztion-1",
    "page": "Pre-Processing",
    "title": "Amplitude Normaliztion",
    "category": "section",
    "text": "Time-domain normalization is used to suppress high-amplitude signals, such as earthquakes or instrumental irregularities. SeisNoise provides time-normalization functions for:one-bit normalization: onebit\nclipping: clip\nrunning mean normalization: smooth\nsuppressing high-amplitude signals: mute"
},

{
    "location": "preprocessing/#Filtering-1",
    "page": "Pre-Processing",
    "title": "Filtering",
    "category": "section",
    "text": "SeisNoise.jl provides bandpass, bandstop, lowpass and highpass filters built from DSP.jl. Due to multiple dispatch in Julia, filter functions in SeisNoise.jl work on either the data in RawData or CorrData objects or directly with Julia Arrays. SeisNoise.jl uses Butterworth filters with a default of 4 corners."
},

{
    "location": "preprocessing/#Filtering-a-RawData-or-CorrData-1",
    "page": "Pre-Processing",
    "title": "Filtering a RawData or CorrData",
    "category": "section",
    "text": "julia> freqmin, freqmax = 1., 10. # low and high frequencies in Hz\njulia> corners = 4 # number of corners in Butterworth filter\njulia> zerophase = true # if true, apply filter twice for no phase shifting\njulia> S = get_data(\"IRIS\",\"TA.V04C..BHZ\",s=\"2006-02-01T00:00:00\",t=\"2006-02-01T01:00:00\")\njulia> SeisIO.demean!(S) # remove mean\njulia> SeisIO.detrend!(S) # remove linear trend\njulia> SeisIO.taper!(S) # taper - defaults to 5% taper on either side of the trace\njulia> bandpass(S,freqmin,freqmax,corners=corners,zerophase=zerophase)"
},

{
    "location": "preprocessing/#Nyquist-Frequency-1",
    "page": "Pre-Processing",
    "title": "Nyquist Frequency",
    "category": "section",
    "text": "Filtering above the Nyquist frequency will give a warning, if using a lowpass or bandpass filter, while using a highpass above the Nyquist frequency will throw an error.julia> S[1].fs / 2. # Nyquist frequency is half sampling rate\n20.\njulia> freqmax = S[1].fs / 2. + 5.\n25.  \njulia> bandpass(S,freqmin,freqmax,corners=corners,zerophase=zerophase)\n┌ Warning: Selected high corner frequency (25.0) of bandpass is at or\n│        above Nyquist (20.0). Applying a high-pass instead.\n└\njulia> lowpass(S,freqmax,corners=corners,zerophase=zerophase)\n┌ Warning: Selected corner frequency (25.0) is\n│ above Nyquist (20.0). Setting Nyquist as high corner.\n└\njulia> highpass(S,freqmax,corners=corners,zerophase=zerophase)\nERROR: frequencies must be less than the Nyquist frequency 20.0"
},

{
    "location": "preprocessing/#SeisIO.demean!-Tuple{AbstractArray{#s38,N} where N where #s38<:AbstractFloat}",
    "page": "Pre-Processing",
    "title": "SeisIO.demean!",
    "category": "method",
    "text": "demean!(A::AbstractArray{<:AbstractFloat})\n\nRemove mean from array A.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisIO.detrend!-Tuple{AbstractArray{#s41,N} where N where #s41<:AbstractFloat}",
    "page": "Pre-Processing",
    "title": "SeisIO.detrend!",
    "category": "method",
    "text": "detrend!(X::AbstractArray{<:AbstractFloat})\n\nRemove linear trend from array X using least-squares regression.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisIO.taper!-Tuple{AbstractArray{#s38,N} where N where #s38<:AbstractFloat,Float64}",
    "page": "Pre-Processing",
    "title": "SeisIO.taper!",
    "category": "method",
    "text": "taper!(A,fs; maxpercentage=0.05, maxlength=20.)\n\nTaper a time series A with samplingrate fs. Defaults to \'hann\' window. Uses smallest of `maxpercentage*fsormax_length`.\n\nArguments\n\nA::AbstractArray: Time series.\nfs::Float64: Sampling rate of time series A.\nmax_percentage::float: Decimal percentage of taper at one end (ranging  from 0. to 0.5).\nmax_length::Float64: Length of taper at one end in seconds.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisNoise.hanningwindow-Tuple{AbstractArray,Int64}",
    "page": "Pre-Processing",
    "title": "SeisNoise.hanningwindow",
    "category": "method",
    "text": "hanningwindow(A,n)\n\nGenerate hanning window of length n.\n\nHanning window is sin(n*pi)^2.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisNoise.bandpass!-Tuple{AbstractArray{#s38,N} where N where #s38<:AbstractFloat,Float64,Float64,Float64}",
    "page": "Pre-Processing",
    "title": "SeisNoise.bandpass!",
    "category": "method",
    "text": "bandpass!(A,freqmin,freqmax,fs,corners=4,zerophase=true)\n\nButterworth-Bandpass Filter.\n\nFilter data A from freqmin to freqmax using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreqmin::Float64: Pass band low corner frequency.\nfreqmax::Float64: Pass band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisNoise.bandstop!-Tuple{AbstractArray{#s38,N} where N where #s38<:AbstractFloat,Float64,Float64,Float64}",
    "page": "Pre-Processing",
    "title": "SeisNoise.bandstop!",
    "category": "method",
    "text": "bandstop!(A,freqmin,freqmax,fs,corners=4,zerophase=true)\n\nButterworth-Bandstop Filter.\n\nFilter data A removing data between frequencies freqmin to freqmax using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreqmin::Float64: Stop band low corner frequency.\nfreqmax::Float64: Stop band high corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisNoise.envelope-Tuple{AbstractArray}",
    "page": "Pre-Processing",
    "title": "SeisNoise.envelope",
    "category": "method",
    "text": "envelope(A)\n\nEnvelope of a function.\n\nComputes the upper and lower envelopes of the given function.\n\nArguments\n\nA::AbstractArray: Data to make envelope of.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisNoise.highpass!-Tuple{AbstractArray{#s38,N} where N where #s38<:AbstractFloat,Float64,Float64}",
    "page": "Pre-Processing",
    "title": "SeisNoise.highpass!",
    "category": "method",
    "text": "highpass(A,freq,fs,corners=4,zerophase=true)\n\nButterworth-Highpass Filter.\n\nFilter data A removing data below certain frequency freq using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisNoise.lowpass!-Tuple{AbstractArray{#s38,N} where N where #s38<:AbstractFloat,Float64,Float64}",
    "page": "Pre-Processing",
    "title": "SeisNoise.lowpass!",
    "category": "method",
    "text": "lowpass(A,freq,fs,corners=4,zerophase=true)\n\nButterworth-Lowpass Filter.\n\nFilter data A over certain frequency freq using corners corners.\n\nArguments\n\nA::AbstractArray: Data to filter\nfreq::Float64: Filter corner frequency.\nfs::Float64: Sampling rate in Hz.\nfs::Int: Filter corners / order.\nzerophase::Bool: If True, apply filter once forwards and once backwards.\n\nThis results in twice the filter order but zero phase shift in the resulting filtered trace.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#SeisNoise.gpufilter!-Tuple{GPUArrays.AbstractGPUArray,DSP.Filters.FilterType,DSP.Filters.ZeroPoleGain}",
    "page": "Pre-Processing",
    "title": "SeisNoise.gpufilter!",
    "category": "method",
    "text": "gpufilter!(A,responsetype,designmethod)\n\nApply filter to array A on the GPU.\n\nArguments\n\nA::AbstractGPUArray: Array on the GPU to filter\nresponsetype::FilterType: DSP.jl filter representation\ndesignmethod::ZeroPoleGain: Filter representation in terms of zeros z, poles p, and\n\ngain k.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#Pre-Processing-Functions-1",
    "page": "Pre-Processing",
    "title": "Pre-Processing Functions",
    "category": "section",
    "text": "Modules = [SeisNoise]\nPages   = [\"ArrayFuncs.jl\", \"filter.jl\"]process_raw!\nonebit!\nclip!\nclamp!\nmute!"
},

{
    "location": "fft/#",
    "page": "Computing FFTs",
    "title": "Computing FFTs",
    "category": "page",
    "text": ""
},

{
    "location": "fft/#Computing-FFTs-1",
    "page": "Computing FFTs",
    "title": "Computing FFTs",
    "category": "section",
    "text": "All correlation in SeisNoise.jl is done in the frequency domain, which can be represented by:C_AB(ω) = u_A(ω) u^_B(ω)Where u_A(ω) is the Fourier transform of ambient noise time series A, u^*_B(ω) is the complex conjugate of the Fourier transform of ambient noise time series B, and C_AB(ω) is the cross-correlation of A and B in the frequency domain. For time and memory efficiency, the real Fourier transform (rfft) is used as opposed to the regular Fourier transform (fft). This gives a speedup of about a factor of 3. A typical workflow for computing fft\'s can include:Real Fourier Transform\nSpectral whitening (removing spectral amplitude information)\nHilbert Transform (for phase-correlation)"
},

{
    "location": "fft/#Creating-FFTData-Structures-1",
    "page": "Computing FFTs",
    "title": "Creating FFTData Structures",
    "category": "section",
    "text": "The rfft function provides the typical workflow for computing ffts in SeisNoise.using SeisNoise, SeisIO\nS = get_data(\"IRIS\",\"TA.V04C..BHZ\",s=\"2006-02-01T00:00:00\",t=\"2006-02-01T01:00:00\")\ncc_step, cc_len = 100., 100.\nR = RawData(S,cc_len,cc_step)\nF = rfft(R)\nFFTData with 35 ffts\n      NAME: \"TA.V04C..BHZ\"                     \n        ID: \"2006-02-01\"                       \n       LOC: 0.0 N, 0.0 E, 0.0 m\n        FS: 40.0\n      GAIN: 1.0\n   FREQMIN: 0.01\n   FREQMAX: 20.0\n    CC_LEN: 100.0\n   CC_STEP: 100.0\n  WHITENED: false                              \n TIME_NORM: \"\"                                 \n      RESP: a0 1.0, f0 1.0, 0z, 0p\n      MISC: 0 entries                          \n     NOTES: 2 entries                          \n         T: 2006-02-01T00:01:40                …\n       FFT: 2001×35 Array{Complex{Float32},2}  Underneath the hood, rfft applies a real Fourier transform to the .x field of a RawData structure, then allocates a new FFTData structure with the Fourier transform data stored in the .fft field."
},

{
    "location": "fft/#Whitening-FFTData-1",
    "page": "Computing FFTs",
    "title": "Whitening FFTData",
    "category": "section",
    "text": "SeisNoise provides three methods for normalizing the spectrum of an FFTData structure. The whiten! method, sets the complex amplitude of frequencies between freqmin and freqmax to 1. This preserves only the phase component of the signal.  freqmin, freqmax = 10., 20.\nwhiten!(F,freqmin,freqmax)The coherence method smooths an the spectrum of an FFTData by the smoothed absolute value of itself.coherence!(F,20)The deconvolution method smooths an the spectrum of an FFTData by the smoothed absolute value squared of itself.deconvolution!(F,20)"
},

{
    "location": "fft/#SeisNoise.whiten!",
    "page": "Computing FFTs",
    "title": "SeisNoise.whiten!",
    "category": "function",
    "text": "whiten!(A, freqmin, freqmax, fs, pad=50)\n\nWhiten spectrum of rfft A between frequencies freqmin and freqmax. Returns the whitened rfft of the time series.\n\nArguments\n\nA::AbstractArray: Time series.\nfs::Real: Sampling rate of time series A.\nfreqmin::Real: Pass band low corner frequency.\nfreqmax::Real: Pass band high corner frequency.\nN::Int: Number of input time domain samples for each rfft.\npad::Int: Number of tapering points outside whitening band.\n\n\n\n\n\nwhiten(F, freqmin, freqmax)\n\nWhiten spectrum of FFTData F between frequencies freqmin and freqmax. Uses real fft to speed up computation. Returns the whitened (single-sided) fft of the time series.\n\nArguments\n\nF::FFTData: FFTData object of fft\'d ambient noise data.\nfreqmin::Real: Pass band low corner frequency.\nfreqmax::Real: Pass band high corner frequency.\npad::Int: Number of tapering points outside whitening band.\n\n\n\n\n\n"
},

{
    "location": "fft/#AbstractFFTs.rfft",
    "page": "Computing FFTs",
    "title": "AbstractFFTs.rfft",
    "category": "function",
    "text": "rfft(R)\n\nComputes windowed rfft of ambient noise data. Returns FFTData structure.\n\nOverloads the rfft function from FFTW.\n\nArguments\n\nR::RawData: RawData structure\n\n\n\n\n\n"
},

{
    "location": "fft/#SeisNoise.coherence!",
    "page": "Computing FFTs",
    "title": "SeisNoise.coherence!",
    "category": "function",
    "text": "coherence!(F,halfwin, waterlevel)\n\nApply coherence method to FFTData F. Where, C_AB(ω) = racu_A(ω) u^_B(ω) u_A(ω)    u_B(ω) \n\nArguments\n\nF::FFTData: FFTData object of fft\'d ambient noise data.\nhalf_win::Int: Number of points in half-window to smooth spectrum.\nwater_level::AbstractFloat: Regularization parameter for spectral smoothing.                               0.01 is a common value [Mehta, 2007].\n\n\n\n\n\n"
},

{
    "location": "fft/#SeisNoise.deconvolution!",
    "page": "Computing FFTs",
    "title": "SeisNoise.deconvolution!",
    "category": "function",
    "text": "deconvolution!(F,halfwin, waterlevel)\n\nApply deconvolution method to FFTData F. Where, C_AB(ω) = racu_A(ω) u^_B(ω) u_B(ω) ^2\n\nArguments\n\nF::FFTData: FFTData object of fft\'d ambient noise data.\nhalf_win::Int: Number of points in half-window to smooth spectrum.\nwater_level::AbstractFloat: Regularization parameter for spectral smoothing.                             0.01 is a common value [Mehta, 2007].\n\n\n\n\n\n"
},

{
    "location": "fft/#SeisNoise.save_fft",
    "page": "Computing FFTs",
    "title": "SeisNoise.save_fft",
    "category": "function",
    "text": "save_fft(F::FFTData, OUT::String)\n\nSave FFTData F to JLD2.\n\n\n\n\n\n"
},

{
    "location": "fft/#SeisNoise.load_fft",
    "page": "Computing FFTs",
    "title": "SeisNoise.load_fft",
    "category": "function",
    "text": "load_fft(filename,chan,day=day)\n\nLoads FFTData for channel chan from JLD2 file filename. If day is specified, loads data from day day, else loads data from all days of chan.\n\n\n\n\n\n"
},

{
    "location": "fft/#Saving/Loading-FFTs-1",
    "page": "Computing FFTs",
    "title": "Saving/Loading FFTs",
    "category": "section",
    "text": "FFTData objects can be saved to disk in the native Julia JLD2 format using the save_fft function.julia> OUTDIR = \"~/TEST/FFT/\"\njulia> save_fft(F,OUTDIR)FFTData are stored in groups by channel (e.g. BHZ or HHZ), then by date (in yyyy-mm-dd format) in JLD2. By default, JLD2 files are saved to /PATH/TO/OUTDIR/NET.STA.CHAN.jld2.file = jldopen(\"~/TEST/FFT/TA.V04C.BHZ.jld2\",\"r\")\nJLDFile ~TEST/FFT/TA.V04C.BHZ.jld2 (read-only)\n └─📂 BHZ\n    └─🔢 2006-02-01To read an FFTData on disk, use the load_fft function:julia> F = load_fft(\"~/TEST/FFT/TA.V04C.BHZ.jld2\",\"BHZ\")\nFFTData with 35 ffts\n      NAME: \"TA.V04C..BHZ\"                     \n        ID: \"2006-02-01\"                       \n       LOC: 0.0 N, 0.0 E, 0.0 m\n        FS: 20.0\n      GAIN: 1.0\n   FREQMIN: 0.05\n   FREQMAX: 5.0\n    CC_LEN: 100                                \n   CC_STEP: 100                                \n  WHITENED: false                              \n TIME_NORM: false                              \n      RESP: c = 0.0, 0 zeros, 0 poles\n      MISC: 0 entries                          \n     NOTES: 2 entries                          \n         T: 2006-02-01T00:00:00.000            …\n       FFT: 2001×35 Array{Complex{Float32},2}Note that it is necessary to specify the channel when using load_fft.whiten!\nSeisNoise.rfft\ncoherence!\ndeconvolution!\nsave_fft\nload_fft"
},

{
    "location": "correlation/#",
    "page": "Correlation",
    "title": "Correlation",
    "category": "page",
    "text": ""
},

{
    "location": "correlation/#Cross-Correlating-1",
    "page": "Correlation",
    "title": "Cross-Correlating",
    "category": "section",
    "text": "Cross-correlation in the frequency domain is element-wise multiplication between u_A(ω), the Fourier transform of ambient noise time series A, and u^*_B(ω), the complex conjugate of the Fourier transform of ambient noise time series B. Options for cross-correlation in SeisNoise.jl includecross-correlation:C_AB(ω) = u_A(ω) u^_B(ω)cross-coherency:C_AB(ω) = fracu_A(ω) u^_B(ω) u_A(omega)    u_B(omega) and deconvolution:C_AB(omega) = fracu_A(omega) u^_B(omega)mid u_B(omega) mid^2The cross-correlation in the time domain is just the inverse real Fourier transform of C_AB(ω):C(τ)_AB = mathfrakF^-1 left(C_AB(ω)right)where τ is the lag time."
},

{
    "location": "correlation/#Computing-Correlations-1",
    "page": "Correlation",
    "title": "Computing Correlations",
    "category": "section",
    "text": "The correlate function provides the typical workflow for computing correlations in SeisNoise.jl. The necessary inputs to correalte are two FFTData structures and the maximum lag time in the correlation in seconds to save, e.g. 200 seconds. Here is an example to cross-correlate data from two Transportable Array stations from February 2, 2006:using SeisNoise, SeisIO\nfs = 40. # sampling frequency in Hz\nfreqmin,freqmax = 0.1,0.2 # minimum and maximum frequencies in Hz\ncc_step, cc_len = 450., 1800. # corrleation step and length in S\nmaxlag = 80. # maximum lag time in correlation\nS1 = get_data(\"IRIS\",\"TA.V04C..BHZ\",s=\"2006-02-01\",t=\"2006-02-02\")\nS2 = get_data(\"IRIS\",\"TA.V05C..BHZ\",s=\"2006-02-01\",t=\"2006-02-02\")\nR1 = RawData(S1,cc_len,cc_step)\nR2 = RawData(S2,cc_len,cc_step)\nF1 = rfft(R1)\nF2 = rfft(R2)\nC = correlate(F1,F2,maxlag)\nCorrData with 188 Corrs\n      NAME: \"TA.V04C..BHZ.TA.V05C..BHZ\"        \n        ID: \"2006-02-01\"                       \n       LOC: 0.0 N, 0.0 E, 0.0 m\n      COMP: \"ZZ\"                               \n   ROTATED: false                              \n CORR_TYPE: \"CC\"                               \n        FS: 40.0\n      GAIN: 1.0\n   FREQMIN: 0.000555556\n   FREQMAX: 20.0\n    CC_LEN: 1800.0\n   CC_STEP: 450.0\n  WHITENED: false                              \n TIME_NORM: \"\"                                 \n      RESP: a0 1.0, f0 1.0, 0z, 0p\n      MISC: 0 entries                          \n     NOTES: 3 entries                          \n      DIST: 0.0\n       AZI: 0.0\n       BAZ: 0.0\n    MAXLAG: 80.0\n         T: 2006-02-01T00:07:30                …\n      CORR: 6401×188 Array{Float32,2}      "
},

{
    "location": "correlation/#SeisNoise.clean_up!",
    "page": "Correlation",
    "title": "SeisNoise.clean_up!",
    "category": "function",
    "text": "clean_up!(A,freqmin,freqmax,fs)\n\nDemean, detrend, taper and filter time series.\n\nArguments\n\nA::AbstractArray: Time series.\nfs::Float64: Sampling rate of time series A in Hz.\nfreqmin::Float64: Pass band low corner frequency in Hz.\nfreqmax::Float64: Pass band high corner frequency in Hz.\n\n\n\n\n\n"
},

{
    "location": "correlation/#SeisNoise.correlate",
    "page": "Correlation",
    "title": "SeisNoise.correlate",
    "category": "function",
    "text": "correlate(FFT1, FFT2, N, maxlag, corr_type=\'cross-correlation\')\n\nCross-correlate ambient noise data in the frequency domain.\n\nCross-correlation can be done using one of three options:\n\nCross-correlation: C_AB(ω) = u_A(ω) u^_B(ω)\nCoherence: C_AB(ω) = racu_A(ω) u^_B(ω) u_A(ω)    u_B(ω) \nDeconvolution: C_AB(ω) = racu_A(ω) u^_B(ω) u_B(ω) ^2\n\nSmoothing of FFTs for coherence and deconvolution should be done before cross-correlating.\n\nArguments\n\nFFT1::AbstractArray: Complex Array of fourier transform of ambient noise data.\nFFT2::AbstractArray: Complex Array of fourier transform of ambient noise data.\nN::Int: Number of input data points in time domain, equal to cc_len * fs.\nmaxlag::Int: Number of data points in cross-correlation to save,                e.g. maxlag = 2000 will save lag times = -2000/fs:2000/fs s.\ncorr_type::String: Type of correlation: cross-correlation, coherence or                      deconv.\n\n\n\n\n\ncorrelate(FFT1, FFT2, maxlag,corr_type=\"CC\")\n\nCross-correlate ambient noise data in the frequency domain.\n\nCross-correlation can be done using one of two options:\n\nCC: Cross-correlation, i.e. C_AB(ω) = u_A(ω) u^_B(ω)\nPCC: Phase cross-correlation, see [Ventosa et al., 2019]\n\nArguments\n\nFFT1::FFTData: FFTData object of fft\'d ambient noise data.\nFFT2::FFTData: FFTData object of fft\'d ambient noise data.\nmaxlag::Float64: Maximum lag time (in seconds) in cross-correlation to save,                    e.g. maxlag = 20. will save lag times = -20.:20. s.\ncorr_type::String: Type of correlation: CC or PCC.\n\n\n\n\n\n"
},

{
    "location": "correlation/#SeisNoise.compute_cc",
    "page": "Correlation",
    "title": "SeisNoise.compute_cc",
    "category": "function",
    "text": "compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;\n           corr_type::String=\"cross-correlation\")\n\nCross-correlate ambient noise data in the frequency domain.\n\nCross-correlation can be done using one of three options:\n\nCross-correlation: C_AB(ω) = u_A(ω) u^_B(ω)\nCoherence: C_AB(ω) = racu_A(ω) u^_B(ω) u_A(ω)    u_B(ω) \nDeconvolution: C_AB(ω) = racu_A(ω) u^_B(ω) u_B(ω) ^2\n\nArguments\n\nFFT1::FFTData: FFTData object of fft\'d ambient noise data.\nFFT2::FFTData: FFTData object of fft\'d ambient noise data.\nmaxlag::Float64: Maximum lag time (in seconds) in cross-correlation to save,                    e.g. maxlag = 20. will save lag times = -20.:20. s.\ncorr_type::String: Type of correlation: cross-correlation, coherence or                      deconv.\n\n\n\n\n\n"
},

{
    "location": "correlation/#SeisNoise.save_corr",
    "page": "Correlation",
    "title": "SeisNoise.save_corr",
    "category": "function",
    "text": "save_corr(C::CorrData, OUT::String)\n\nSave CorrData C to JLD2.\n\n\n\n\n\n"
},

{
    "location": "correlation/#SeisNoise.load_corr",
    "page": "Correlation",
    "title": "SeisNoise.load_corr",
    "category": "function",
    "text": "load_corr(filename,chan,day=day)\n\nLoads CorrData for channel chan on day day from JLD2 file filename.\n\n\n\n\n\n"
},

{
    "location": "correlation/#SeisNoise.stack!",
    "page": "Correlation",
    "title": "SeisNoise.stack!",
    "category": "function",
    "text": "stack!(C)\n\nStack correlation by time interval. The default is to stack by day. Using allstack == true will stack all available correlations. To use phase-weighted stack, specify the amount of phase_smoothing in seconds.\n\nArguments\n\nC::CorrData: Correlation data.\ninterval::Period: Interval over which to stack C.\nallstack::Bool: If true, stack all data.\nstack_type::Function: Type of stacking. Options are mean, pws, robuststack, etc..\n\n\n\n\n\n"
},

{
    "location": "correlation/#Saving/Loading-Correlations-1",
    "page": "Correlation",
    "title": "Saving/Loading Correlations",
    "category": "section",
    "text": "CorrData objects can be saved to disk in the native Julia JLD2 format using the save_corr function.julia> OUTDIR = \"~/TEST/CORR/\"\njulia> save_corr(C,OUTDIR)CorrData are stored in groups by component (e.g. ZZ or RZ), then by date (in yyyy-mm-dd format) in JLD2. By default, JLD2 files are saved to /PATH/TO/OUTDIR/NET1.STA1.CHAN1.NET2.STA2.CHAN2.jld2.file = jldopen(\"~/TEST/CORR/TA.V04C.BHZ.TA.V05C.BHZ.jld2\",\"r\")\nJLDFile ~/TEST/CORR/TA.V04C.BHZ.TA.V05C.BHZ.jld2 (read-only)\n └─📂 ZZ\n    └─🔢 2006-02-01To read an CorrData on disk, use the load_corr function:julia> C = load_corr(\"~/TEST/CORR/TA.V04C.BHZ.TA.V05C.BHZ.jld2\",\"ZZ\")\nCorrData with 188 Corrs\n      NAME: \"TA.V04C..BHZ.TA.V05C..BHZ\"        \n        ID: \"2006-02-01\"                       \n       LOC: 0.0 N, 0.0 E, 0.0 m\n      COMP: \"ZZ\"                               \n   ROTATED: false                              \n CORR_TYPE: \"coherence\"                        \n        FS: 40.0\n      GAIN: 1.0\n   FREQMIN: 0.1\n   FREQMAX: 0.2\n    CC_LEN: 1800                               \n   CC_STEP: 450                                \n  WHITENED: false                              \n TIME_NORM: false                              \n      RESP: c = 0.0, 0 zeros, 0 poles\n      MISC: 0 entries                          \n     NOTES: 2 entries                          \n    MAXLAG: 80.0\n         T: 2006-02-01T00:07:30.000            …\n      CORR: 6401×188 Array{Float32,2}  clean_up!\ncorrelate\ncompute_cc\nsave_corr\nload_corr\nstack!"
},

{
    "location": "extend/#",
    "page": "Extending SeisNoise",
    "title": "Extending SeisNoise",
    "category": "page",
    "text": ""
},

{
    "location": "extend/#Extending-SeisNoise-1",
    "page": "Extending SeisNoise",
    "title": "Extending SeisNoise",
    "category": "section",
    "text": "Extending SeisNoise for your ambient noise workflow should be fairly easy. Let\'s say you have a function written in Julia that you want to apply to a RawData structure. For example, below is a function double! that accepts an array of Real numbers and doubles them in-place.function double!(A::AbstractArray{T}) where T <: Real\n    A .*= T(2)\n    return nothing\nendIf you want to apply double! to the data in a RawData or CorrData structure, you need to input R.x or C.corr to the double! function, as shown below:using SeisNoise\nR = RawData()\nA = rand(Float32,6,4)\nR.x = deepcopy(A)\ndouble!(R.x)\nR.x\n6×4 Array{Float32,2}:\n 1.78214    1.02787   0.55617   0.512943\n 1.35262    1.54597   0.212465  1.7978\n 1.94816    1.53011   1.6552    0.795328\n 0.955228   0.787446  1.43849   0.175546\n 0.714791   1.05514   0.124099  1.51923\n 0.0911338  1.22735   1.55351   1.15741Rather than inputting R.x, one could use multiple-dispatch to define another double! function that accepts RawData, like this:double!(R::RawData) = double!(R.x)Now we can input a RawData structure to our double! function, like soR.x = deepcopy(A)\ndouble!(R)\nR.x\n6×4 Array{Float32,2}:\n 1.78214    1.02787   0.55617   0.512943\n 1.35262    1.54597   0.212465  1.7978\n 1.94816    1.53011   1.6552    0.795328\n 0.955228   0.787446  1.43849   0.175546\n 0.714791   1.05514   0.124099  1.51923\n 0.0911338  1.22735   1.55351   1.15741By convention, the ! in double! implies that the operation is applied in-place, meaning that the input array is overwritten and no new memory is allocated. Allocating versions of double! could look like this:double(A::AbstractArray) = (U = deepcopy(A); double!(U); return U)\ndouble(R::RawData) = (U = deepcopy(R); double!(U); return U)So now we can output a new RawData structure with the allocating double:R.x = deepcopy(A)\nRnew = double(R)\nRnew.x\n6×4 Array{Float32,2}:\n 1.78214    1.02787   0.55617   0.512943\n 1.35262    1.54597   0.212465  1.7978\n 1.94816    1.53011   1.6552    0.795328\n 0.955228   0.787446  1.43849   0.175546\n 0.714791   1.05514   0.124099  1.51923\n 0.0911338  1.22735   1.55351   1.15741"
},

{
    "location": "extend/#Developer-Advice-1",
    "page": "Extending SeisNoise",
    "title": "Developer Advice",
    "category": "section",
    "text": "We recommend first writing kernel functions that accept AbstractArrays, as shown with double!, then writing a wrapper function that accepts RawData, FFTData, or CorrData objects. Writing code in this way leads to 1) faster code 2) code reuse and 3) type-stability. If you are interested in writing high-performance code, we recommend having a look at the Julia Performance Tips."
},

{
    "location": "extend/#Adding-Methods-to-SeisNoise-1",
    "page": "Extending SeisNoise",
    "title": "Adding Methods to SeisNoise",
    "category": "section",
    "text": "If you have a method/routine for processing ambient noise cross-correlations that you think would be helpful for the community, feel free to let us know. Check out the Contributor\'s Guide to learn how to add your code/ideas to SeisNoise. "
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "Coming soon to a doc page near you!"
},

{
    "location": "contributing/#",
    "page": "Contributer\'s Guide",
    "title": "Contributer\'s Guide",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/#Contributor\'s-Guide-1",
    "page": "Contributer\'s Guide",
    "title": "Contributor\'s Guide",
    "category": "section",
    "text": "Thank you for considering contributing to SeisNoise! This short guide will give you ideas on how you can contribute and help you make a contribution.Please feel free to ask us questions and chat with us at any time if you\'re unsure about anything."
},

{
    "location": "contributing/#What-can-I-do?-1",
    "page": "Contributer\'s Guide",
    "title": "What can I do?",
    "category": "section",
    "text": "Try to using SeisNoise process one or two days of cross-correlations. If you run into any problems or find it difficult to use or understand, please open an issue!\nWrite up an example or tutorial on how to do something useful with SeisNoise, like processing cross-correlation results on the GPU.\nImprove documentation or comments if you found something hard to use.\nImplement a new feature if you need it to use SeisNoise.If you\'re interested in working on something, let us know by commenting on existing issues or by opening a new issue if. This is to make sure no one else is working on the same issue and so we can help and guide you in case there is anything you need to know beforehand."
},

{
    "location": "contributing/#Ground-Rules-1",
    "page": "Contributer\'s Guide",
    "title": "Ground Rules",
    "category": "section",
    "text": "Each pull request should consist of a logical collection of changes. You can include multiple bug fixes in a single pull request, but they should be related. For unrelated changes, please submit multiple pull requests.\nDo not commit changes to files that are irrelevant to your feature or bugfix (eg: .gitignore).\nBe willing to accept criticism and work on improving your code; we don\'t want to break other users\' code, so care must be taken not to introduce bugs. We discuss pull requests and keep working on them until we believe we\'ve done a good job.\nBe aware that the pull request review process is not immediate, and is generally proportional to the size of the pull request."
},

{
    "location": "contributing/#Reporting-a-bug-1",
    "page": "Contributer\'s Guide",
    "title": "Reporting a bug",
    "category": "section",
    "text": "The easiest way to get involved is to report issues you encounter when using SeisNoise or by requesting something you think is missing.Head over to the issues page.\nSearch to see if your issue already exists or has even been solved previously.\nIf you indeed have a new issue or request, click the \"New Issue\" button.\nPlease be as specific as possible. Include the version of the code you were using, as well as what operating system you are running. The output of Julia\'s versioninfo() and ] status is helpful to include. If possible, include complete, minimal example code that reproduces the problem."
},

{
    "location": "contributing/#Setting-up-your-development-environment-1",
    "page": "Contributer\'s Guide",
    "title": "Setting up your development environment",
    "category": "section",
    "text": "Install Julia on your system.\nInstall git on your system if it is not already there (install XCode command line tools on a Mac or git bash on Windows).\nLogin to your GitHub account and make a fork of the SeisNoise repository by clicking the \"Fork\" button.\nClone your fork of the SeisNoise repository (in terminal on Mac/Linux or git shell/ GUI on Windows) in the location you\'d like to keep it.\ngit clone https://github.com/your-user-name/SeisNoise.jl.git\nNavigate to that folder in the terminal or in Anaconda Prompt if you\'re on Windows.\nConnect your repository to the upstream (main project).\ngit remote add SeisNoise https://github.com/tclements/SeisNoise.jl.git\nCreate the development environment by opening Julia via julia --project then typing in ] instantiate. This will install all the dependencies in the Project.toml file.\nYou can test to make sure SeisNoise works by typing in ] test which will run all the tests (this can take a while).Your development environment is now ready!"
},

{
    "location": "contributing/#Pull-Requests-1",
    "page": "Contributer\'s Guide",
    "title": "Pull Requests",
    "category": "section",
    "text": "Changes and contributions should be made via GitHub pull requests against the master branch.When you\'re done making changes, commit the changes you made. Chris Beams has written a guide on how to write good commit messages.When you think your changes are ready to be merged into the main repository, push to your fork and submit a pull request.Working on your first Pull Request? You can learn how from this free video series How to Contribute to an Open Source Project on GitHub, Aaron Meurer\'s tutorial on the git workflow, or the guide “How to Contribute to Open Source\"."
},

{
    "location": "contributing/#Documentation-1",
    "page": "Contributer\'s Guide",
    "title": "Documentation",
    "category": "section",
    "text": "Now that you\'ve made your awesome contribution, it\'s time to tell the world how to use it. Writing documentation strings is really important to make sure others use your functionality properly. Didn\'t write new functions? That\'s fine, but be sure that the documentation for the code you touched is still in great shape. It is not uncommon to find some strange wording or clarification that you can take care of while you are here."
},

{
    "location": "contributing/#Credits-1",
    "page": "Contributer\'s Guide",
    "title": "Credits",
    "category": "section",
    "text": "This contributor\'s guide is heavily based on the Oceananigans contributor\'s guide, which is based upon the MetPy contributor\'s guide."
},

{
    "location": "func_index/#",
    "page": "Function Index",
    "title": "Function Index",
    "category": "page",
    "text": ""
},

{
    "location": "func_index/#Index-1",
    "page": "Function Index",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
