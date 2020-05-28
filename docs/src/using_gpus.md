# Using GPUs

SeisNoise is written to run on graphical processing units (GPUs)
for increased performance. Depending on your CPU and GPU combination, speedups of >10-20x
are possible.

!!! tip "Running on GPUs"
    If you are having issues with running SeisNoise on a GPU or setting things up,
    please [open an issue](https://github.com/tclements/SeisNoise.jl/issues/new)
    and we'll do our best to help out!

## When to use a GPU

GPUs are very useful for number crunching. If you have a large number of stations to cross-correlate, using a GPU *could* greatly improve your cross-correlation processing time.

GPU processing tends to be memory-limited. High-end GPUs (think $$$$) such as the
Nvidia Tesla V100 only come with up to 32 GB of memory. On a GPU with 16 GB of memory, you can cross-correlate ~100 day-long 3-component stations sampled at 100 Hz in memory.

## Getting access to GPUs

If you don't have a GPU there are a few resources you can try to acquire a GPU from.

If you have access to any supercomputer clusters, check to see if they have any GPUs.
See also this Stack Overflow post:
[Where can I get access to GPU cluster for educational purpose?](https://scicomp.stackexchange.com/questions/8508/where-can-i-get-access-to-gpu-cluster-for-educational-purpose)

Cloud computing providers such as Google Cloud and Amazon EC2 allow you to rent GPUs per
hour. Sometimes they offer free trials or credits that can be used towards GPUs although
they seem to be getting less common.

See the [Julia on Google Colab: Free GPU-Accelerated Shareable Notebooks](https://discourse.julialang.org/t/julia-on-google-colab-free-gpu-accelerated-shareable-notebooks/15319)
post on the Julia Discourse.

[Code Ocean](https://codeocean.com/) also has
[GPU support](https://help.codeocean.com/en/articles/1053107-gpu-support) and allows you
to spin up capsules with pretty decent Tesla K80 GPUs for free (for now) if you want to
play around with them.  You'll
want to use their "Ubuntu Linux with GPU support (18.04.3)" with the ability to compile
CUDA code. Then you'll have to install Julia manually.

[NextJournal](https://nextjournal.com/) provides a Jupyter-like Julia inferface with access to Tesla K80 GPUs for free.

## I have a GPU. Now what?

Make sure you have an Nvidia GPU that is CUDA compatible:
[https://developer.nvidia.com/cuda-gpus](https://developer.nvidia.com/cuda-gpus). Most
recent GPUs should be but older GPUs and many laptop GPUs may not be.

Then download and install the CUDA Toolkit:
[https://developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads)

Once the CUDA Toolkit is installed, you might have to build SeisNoise again
```
julia>]
(v1.4) pkg> build SeisNoise
```

## Using the GPU with Julia

SeisNoise uses the [Cuda.jl library](https://github.com/JuliaGPU/CUDA.jl) to do all GPU-based operations. To get started with GPU-programming, take a look at the [Cuda.jl introduction tutorial](https://juliagpu.gitlab.io/CUDA.jl/tutorials/introduction/). Here is a great slide from Tim Besard on [GPU computing in Julia](https://docs.google.com/presentation/d/1l-BuAtyKgoVYakJSijaSqaTL3friESDyTOnU2OLqGoA/edit#slide=id.p).

Allocating arrays on the GPU is simple using the `cu` function. This example creates a random maxtrix on the CPU, then sends it to the GPU using `cu`:

```julia
using CUDA
A = cu(rand(1000,10))
1000×10 CuArray{Float32,2,Nothing}:
 0.483598   0.0278827  0.564381  …  0.0179997  0.0894625  0.943952
 0.532773   0.601974   0.965769     0.91777    0.746491   0.837131
 0.0613844  0.638339   0.401549     0.385933   0.627327   0.682377
 0.597944   0.349372   0.147358     0.607321   0.592268   0.381262
 0.516212   0.505277   0.35825      0.0874267  0.946768   0.525376
 ⋮                               ⋱                        
 0.696375   0.694293   0.646136     0.0734224  0.0999966  0.804346
 0.266413   0.710535   0.762904     0.13076    0.0786623  0.914678
 0.255818   0.0446712  0.943867     0.981916   0.835427   0.506067
 0.611743   0.715961   0.181798     0.793934   0.568452   0.84297

```

or similarly using the Julia pipe syntax,

```julia
A = rand(1000,10) |> cu
1000×10 CuArray{Float32,2,Nothing}:
 0.511658   0.231007   0.819      …  0.110045   0.517574    0.195391
 0.589017   0.69037    0.0288187     0.599881   0.0290229   0.136348
 0.175985   0.217625   0.124626      0.0619627  0.862026    0.431731
 0.0430858  0.569207   0.585091      0.605072   0.333407    0.505581
 0.510674   0.969721   0.875608      0.984235   0.0604092   0.0494258
 ⋮                                ⋱                         
 0.938167   0.673817   0.0642641     0.401075   0.420085    0.0400348
 0.311241   0.0985884  0.309681      0.33953    0.793376    0.385521
 0.738205   0.500696   0.294173      0.681627   0.00480331  0.751362
 0.740553   0.566934   0.845259      0.503732   0.970519    0.218747

```
Note that here `A` is a `CuArray` which is stored in memory on the GPU.
Most standard Julia functions (matrix multiplication, FFT, mean, max, etc..) will work out of the box on `CuArray`s. For example, here is matrix multiplication

```julia
B = rand(10,1000) |> cu
A * B
1000×1000 CuArray{Float32,2,Nothing}:
 3.02011  2.85182  3.17944  …  2.20066  1.89412  2.32329
 3.12874  2.83144  3.34268     1.65363  1.97212  2.41745
 3.22296  2.53652  2.84638     2.11542  2.31627  1.71027
 3.21165  2.51493  2.9726      2.07334  2.55837  2.63839
 3.2731   2.83646  3.75265     1.44419  1.68725  3.5038
 2.18507  1.80968  2.13033  …  1.47744  1.97971  1.59087
 ⋮                          ⋱                    
 1.95406  1.95351  2.62476  …  1.07703  1.49473  1.72594
 2.26404  2.09106  2.86892     1.12002  1.21775  1.73231
 2.56381  2.1512   2.70032     2.03604  1.98124  1.38784
 3.07007  3.02482  3.36264     2.14246  2.4654   2.34386
 3.46125  2.72799  3.75003     1.8953   2.44945  2.65501
```

## SeisNoise on the GPU
SeisNoise can process data and compute cross-correlations on the GPU with CUDA. CUDA.jl provides an the (`CuArray`) type for storing data on the GPU. Data in SeisNoise structures (`R.x`, `F.fft`, and `C.corr` fields, for `RawData`, `FFTData`, and `CorrData`, respectively) can move between an `Array` on the CPU to a `CuArray` on the GPU using the `gpu` and `cpu` functions, as shown below.   

!!! note
    Only **Nvidia** GPUs are suported at the moment. Hold in there for AMD/OpenCL/Metal support...

```julia
# create raw data and send to GPU
R = RawData(S1, cc_len, cc_step) |> gpu
R.x
72000×188 CuArrays.CuArray{Float32,2,Nothing}

# send data back to the CPU
R = R |> cpu
R.x
72000×188 Array{Float32,2}
```

All basic processing remains the same on the GPU. Here is a complete cross-correlation routine on the GPU:

```julia
# send data to GPU
R1 = RawData(S1, cc_len, cc_step) |> gpu
R2 = RawData(S2, cc_len, cc_step) |> gpu
R = [R1,R2]

# preprocess on the GPU
detrend!.(R)
taper!.(R)
bandpass!.(R,freqmin,freqmax,zerophase=true)

# Real FFT on GPU
FFT = rfft.(R)
whiten!.(FFT,freqmin,freqmax)

# compute correlation and send to cpu
C = correlate(FFT[1],FFT[2],maxlag) |> cpu
```

### Routines Implemented on the GPU
- Preprocessing:
  - `detrend`,`demean`, `taper`, `onebit`, `smooth`
- Filtering:
  - `bandpass`, `bandstop`, `lowpass`, `highpass`
- Fourier Domain:
  - `whiten`, `rfft`, `irfft`
- Cross-correlation:
  - `correlate`, `cross-coherence`, `deconvolution`
- Post-processing:
  - `stack`, `filter`s, etc..


### Credits

This page builds heavily on the using GPUs page from [Oceananigans](https://clima.github.io/Oceananigans.jl/stable/using_gpus/).
