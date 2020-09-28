import SeisIO: demean, demean!, taper, taper!,detrend, detrend!
export detrend, detrend!, taper, taper!, demean, demean!, hanningwindow
# Signal processing functions for arrays (rather than SeisData or SeisChannel)

"""
    detrend!(X::AbstractArray{<:AbstractFloat})

Remove linear trend from array `X` using least-squares regression.
"""
function detrend!(X::AbstractArray{<:AbstractFloat})
    T = eltype(X)
    N = size(X,1)

    # create linear trend matrix
    A = similar(X,T,N,2)
    A[:,2] .= T(1)
    A[:,1] .= range(T(0),T(1),length=N)
    # create linear trend matrix
    R = transpose(A) * A

    # do the matrix inverse for 2x2 matrix
    # this is really slow on GPU
    Rinv = inv(Array(R)) |> typeof(R)
    factor = Rinv * transpose(A)

    # remove trend
    X .-= A * (factor * X)
    return nothing
end
detrend(A::AbstractArray{<:AbstractFloat}) = (U = deepcopy(A);
        detrend!(U);return U)
detrend!(R::RawData) = detrend!(R.x)
detrend(R::RawData) = (U = deepcopy(R); detrend!(U.x); return U)
detrend!(C::CorrData) = detrend!(C.corr)
detrend(C::CorrData) = (U = deepcopy(C); detrend!(U.corr); return U)


"""
    demean!(A::AbstractArray{<:AbstractFloat})

Remove mean from array `A`.
"""
function demean!(A::AbstractArray{<:AbstractFloat}; dims=1)
      μ = mean(A,dims=dims)
      A .-= μ
  return nothing
end
demean(A::AbstractArray{<:AbstractFloat}; dims=1) = (U = deepcopy(A);
       demean!(U,dims=dims);return U)
demean!(R::RawData) = demean!(R.x)
demean(R::RawData) = (U = deepcopy(R); demean!(U.x); return U)
demean!(C::CorrData) = demean!(C.corr)
demean(C::CorrData) = (U = deepcopy(C); demean!(U.corr); return U)

"""
   taper!(A,fs; max_percentage=0.05, max_length=20.)

Taper a time series `A` with sampling_rate `fs`.
Defaults to 'hann' window. Uses smallest of `max_percentage` * `fs`
or `max_length`.

# Arguments
- `A::AbstractArray`: Time series.
- `fs::AbstractFloat`: Sampling rate of time series `A`.
- `max_percentage::float`: Decimal percentage of taper at one end (ranging
   from 0. to 0.5).
- `max_length::Real`: Length of taper at one end in seconds.
"""
function taper!(A::AbstractArray{<:AbstractFloat}, fs::Real;
                max_percentage::AbstractFloat=0.05, max_length::Real=20.)
   Nrows = size(A,1)
   T = eltype(A)
   wlen = min(Int(floor(Nrows * max_percentage)), Int(floor(max_length * fs)), Int(
         floor(Nrows/2)))
   taper_sides = hanningwindow(A,2 * wlen)
   A[1:wlen,:] .*= taper_sides[1:wlen]
   A[end-wlen:end,:] .*= taper_sides[wlen:end]
   return nothing
end
taper(A::AbstractArray{<:AbstractFloat}, fs::Real;
       max_percentage::AbstractFloat=0.05,max_length::Real=20.) = (U = deepcopy(A);
       taper!(U,fs,max_percentage=max_percentage,max_length=max_length);return U)
taper!(R::RawData; max_percentage::AbstractFloat=0.05,
       max_length::Real=20.) = taper!(R.x,R.fs,max_percentage=max_percentage,
       max_length=max_length)
taper(R::RawData; max_percentage::AbstractFloat=0.05,
       max_length::Real=20.) = (U = deepcopy(R); taper!(U.x,U.fs,
       max_percentage=max_percentage,max_length=max_length); return U)
taper!(C::CorrData; max_percentage::AbstractFloat=0.05,
      max_length::Real=20.) = taper!(C.corr,C.fs,max_percentage=max_percentage,
      max_length=max_length)
taper(C::CorrData; max_percentage::AbstractFloat=0.05,
      max_length::Real=20.) = (U = deepcopy(C); taper!(U.corr,U.fs,
      max_percentage=max_percentage,max_length=max_length); return U)


"""
  hanningwindow(A,n)

Generate hanning window of length `n`.

Hanning window is sin(n*pi)^2.
"""
function hanningwindow(A::AbstractArray, n::Int)
   T = eltype(A)
   win = similar(A,T,n)
   win .= T(pi) .* range(T(0.),stop=T(n),length=n)
   win .= sin.(win).^2
   return win
end
