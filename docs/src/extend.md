# Extending SeisNoise

Extending SeisNoise for your ambient noise workflow should be fairly easy. Let's say you
have a function written in Julia that you want to apply to a `RawData` structure.
For example, below is a function `double!` that accepts an array of `Real` numbers
and doubles them in-place.


```julia
function double!(A::AbstractArray{T}) where T <: Real
    A .*= T(2)
    return nothing
end
```

If you want to apply `double!` to the data in a `RawData` or `CorrData` structure, you need
to input `R.x` or `C.corr` to the `double!` function, as shown below:

```julia
using SeisNoise
R = RawData()
A = rand(Float32,6,4)
R.x = deepcopy(A)
double!(R.x)
R.x
6×4 Array{Float32,2}:
 1.78214    1.02787   0.55617   0.512943
 1.35262    1.54597   0.212465  1.7978
 1.94816    1.53011   1.6552    0.795328
 0.955228   0.787446  1.43849   0.175546
 0.714791   1.05514   0.124099  1.51923
 0.0911338  1.22735   1.55351   1.15741
```

Rather than inputting `R.x`, one could use [multiple-dispatch](https://en.wikipedia.org/wiki/Multiple_dispatch)
to define another `double!` function that accepts `RawData`, like this:

```julia
double!(R::RawData) = double!(R.x)
```

Now we can input a `RawData` structure to our `double!` function, like so

```julia
R.x = deepcopy(A)
double!(R)
R.x
6×4 Array{Float32,2}:
 1.78214    1.02787   0.55617   0.512943
 1.35262    1.54597   0.212465  1.7978
 1.94816    1.53011   1.6552    0.795328
 0.955228   0.787446  1.43849   0.175546
 0.714791   1.05514   0.124099  1.51923
 0.0911338  1.22735   1.55351   1.15741
```

By convention, the `!` in `double!` implies that the operation is applied in-place,
meaning that the input array is overwritten and no new memory is allocated.
Allocating versions of `double!` could look like this:

```julia
double(A::AbstractArray) = (U = deepcopy(A); double!(U); return U)
double(R::RawData) = (U = deepcopy(R); double!(U); return U)
```

So now we can output a new `RawData` structure with the allocating `double`:

```julia
R.x = deepcopy(A)
Rnew = double(R)
Rnew.x
6×4 Array{Float32,2}:
 1.78214    1.02787   0.55617   0.512943
 1.35262    1.54597   0.212465  1.7978
 1.94816    1.53011   1.6552    0.795328
 0.955228   0.787446  1.43849   0.175546
 0.714791   1.05514   0.124099  1.51923
 0.0911338  1.22735   1.55351   1.15741
```

## Developer Advice
We recommend first writing kernel functions that accept `AbstractArray`s, as shown
with `double!`, then writing a wrapper function that accepts `RawData`, `FFTData`,
or `CorrData` objects. Writing code in this way leads to 1) faster code 2) code reuse
and 3) type-stability. If you are interested in writing high-performance code,
we recommend having a look at the [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/index.html).

## Adding Methods to SeisNoise
If you have a method/routine for processing ambient noise cross-correlations that
you think would be helpful for the community, feel free to let us know. Check out the
[Contributor's Guide](@ref) to learn how to add your code/ideas to SeisNoise. 
