export stack, stack!, remove_nan, remove_nan!, shorten, shorten!
export pws, pws!, robuststack, robuststack!, adaptive_filter!, adaptive_filter
export robustpws, robustpws!
"""

  stack!(C)

Stack correlation by time interval.
The default is to stack by day. Using `allstack == true` will stack all available
correlations. To use phase-weighted stack, specify the amount of
`phase_smoothing` in seconds.

# Arguments
- `C::CorrData`: Correlation data.
- `interval::Period`: Interval over which to stack `C`.
- `allstack::Bool`: If `true`, stack all data.
- `stack_type::Function`: Type of stacking. Options are mean, pws, robuststack, etc..
"""
function stack!(C::CorrData; interval::Period=Day(1),
                allstack::Bool=false,stacktype::Function=mean)

    if allstack == true
          if stacktype == mean
                C.corr = mean(C.corr,dims=2)
          else # use phase-weighted stack
                C.corr = stacktype(C.corr)
          end
          C.t = [C.t[1]]
    else # stack by interval
          C.t = d2u.(round.(u2d.(C.t),interval,RoundDown))
          stackT = unique(C.t)
          ind = indexin(stackT,C.t)
          push!(ind,length(C.t)+1)
          stack_out = Array{eltype(C.corr)}(undef,size(C.corr,1),length(stackT))

          for ii = 1:length(stackT)
                 if stacktype == mean
                       stack_out[:,ii] = mean(C.corr[:,ind[ii]:ind[ii+1]-1],dims=2)
                 else
                       stack_out[:,ii] = stacktype(C.corr[:,ind[ii]:ind[ii+1]-1])
                 end
          end

          C.t = stackT
          C.corr = stack_out
    end
    return C
end
stack(C::CorrData; interval::Period=Day(1),
      allstack::Bool=false,stacktype::Function=mean) = (U = deepcopy(C);
      stack!(U,interval=interval,allstack=allstack,
             stacktype=stacktype);return U)

"""

  robuststack(A)

Performs robust stack on array `A`.

Follows methods of Pavlis and Vernon, 2010.

# Arguments
- `A::AbstractArray`: Time series stored in columns.
- `ϵ::AbstractFloat`: Threshold for convergence of robust stack.
- `maxiter::Int`: Maximum number of iterations to converge to robust stack.
"""
function robuststack(A::AbstractArray{T};ϵ::AbstractFloat=Float32(1e-4),
                     maxiter::Int=10) where T <: AbstractFloat
    N = size(A,2)
    Bold = median(A,dims=2)
    Bold ./= norm(Bold,2)
    w = Array{T}(undef,N)
    r = Array{T}(undef,N)
    d2 = Array{T}(undef,N)

    # do 2-norm for all columns in A
    for ii = 1:N
        d2[ii] = norm(A[:,ii],2)
    end

    BdotD = sum(A .* Bold,dims=1)

    for ii = 1:N
        r[ii] = norm(A[:,ii] .- (BdotD[ii] .* Bold),2)
        w[ii] = abs(BdotD[ii]) ./ d2[ii] ./ r[ii]
    end

    Bnew = mean(A,weights(w),dims=2)
    Bnew ./= norm(Bnew,2)

    # check convergence
    ϵN = norm(Bnew .- Bold,1) / (norm(Bnew,2) * N)
    Bold = Bnew
    iter = 0
    while (ϵN > ϵ) && (iter <= maxiter)
        BdotD = sum(A .* Bold,dims=1)

        for ii = 1:N
            r[ii] = norm(A[:,ii] .- (BdotD[ii] .* Bold),2)
            w[ii] = abs(BdotD[ii]) ./ d2[ii] ./ r[ii]
        end

        Bnew = mean(A,weights(w),dims=2)
        Bnew ./= norm(Bnew,2)

        # check convergence
        ϵN = norm(Bnew .- Bold,1) / (norm(Bnew,2) * N)
        Bold = Bnew
        iter += 1
    end
    return Bnew
end
robuststack!(C::CorrData;ϵ::AbstractFloat=eps(Float32),maxiter::Int=10) =
       (C.corr = robuststack(C.corr,ϵ=ϵ,maxiter=maxiter); C.t = C.t[1:1]; return C)
robuststack(C::CorrData;ϵ::AbstractFloat=eps(Float32),maxiter::Int=10)  =
       (U = deepcopy(C); U.corr = robuststack(U.corr,ϵ=ϵ,maxiter=maxiter); U.t = U.t[1:1];
       return U)

"""

  pws(A)

Performs phase-weighted stack on array `A` of time series.

Follows methods of Schimmel and Paulssen, 1997.
If s(t) is time series data,
S(t) = s(t) + i*H(s(t)), where H(s(t)) is Hilbert transform of s(t)
S(t) = s(t) + i*H(s(t)) = A(t)*exp(i*phi(t)), where
A(t) is envelope of s(t) and phi(t) is phase of s(t)
Phase-weighted stack, g(t), is then:
g(t) = 1/N sum j = 1:N s_j(t) * | 1/N sum k = 1:N exp[i * phi_k(t)]|^v
where N is number of traces used, v is sharpness of phase-weighted stack

# Arguments
- `A::AbstractArray`: Time series stored in columns.
- `pow::Int`: Sharpness of transition from phase similarity to dissimilarity.
"""
function pws(A::AbstractArray{T}; pow::Real=2) where T <: AbstractFloat
    # preserve type-stability
    if !isa(pow,Int)
        pow = T(pow)
    end
    Nrows,Ncols = size(A)
    phase_stack = abs.(sum(exp.(im .* angle.(hilbert(A))),dims=2) ./ Ncols) .^ pow
    return mean(A .* phase_stack,dims=2)
end
pws!(C::CorrData; pow::Real=2) = (C.corr=pws(C.corr,pow=pow); C.t = C.t[1:1]; return nothing)
pws(C::CorrData; pow::Real=2) = (U = deepcopy(C);
    U.corr=pws(U.corr,pow=pow); U.t = U.t[1:1];return U)

"""

  remove_nan!(C)

Remove correlations with Nan from CorrData `C`.
"""
function remove_nan!(C::CorrData)
      ind = []
      for ii = 1:length(C.t)
            if !any(isnan,C.corr[:,ii])
                  append!(ind,ii)
            end
      end

      # throw error if all values are NaN
      if length(ind) == 0
          throw(ArgumentError("All correlation windows contain NaNs."))
      end

      # return non-NaN columns
      C.corr = C.corr[:,ind]
      C.t = C.t[ind]
      return nothing
end
remove_nan(C::CorrData) = (U = deepcopy(C);remove_nan!(U);return U)

"""

  shorten!(C)

Clip CorrData `C` from lags τ = 0 to abs(maxlag).
"""
function shorten!(C::CorrData, maxlag::Real)
      if maxlag >= C.maxlag || maxlag <= 0.
            return C
      end

      # get timearray
      lags = -C.maxlag:1/C.fs:C.maxlag
      ind = findall(x -> abs(x) <= maxlag, lags)
      C.maxlag = maxlag
      C.corr = C.corr[ind,:]
      return nothing
end
shorten(C::CorrData,maxlag::Real) = (U = deepcopy(C); shorten!(U,maxlag);return U)

"""

  adaptive_filter!(A,window,fs)

Adaptive covariance filter to enhance coherent signals.

Fellows the method of Nakata et al., 2015 (Appendix B). The filtered signal
`X` is given by `X = irfft(P*X(w))`` where `X` is the FFT'ed spectra of `A`
and `P` is the filter. `P` is constructed by using the temporal covariance matrix.
​
# Arguments
- `A::AbstractArray`: Array of daily/hourly cross-correlation functions
- `window::AbstractFloat`: Window length to apply adaptive filter [s].
- `fs::Real`: Sampling rate of time series `A`.
- `g::Int`: Positive number to adjust the filter harshness
- `overlap::AbstractFloat`: Percent overlap between windows
"""
function adaptive_filter!(A::AbstractArray{T}, window::AbstractFloat,
                          fs::Real; g::AbstractFloat=2., overlap::AbstractFloat=0.9) where T <: AbstractFloat
    if ndims(A) == 1
        return nothing
    end

    if overlap <= 0 || overlap >= 1
        throw(ArgumentError("overlap must be > 0 and < 1."))
    end

    Nrows, Ncols = size(A)
    window_samples = convert(Int,round(window * fs,digits=0))
    window_step = convert(Int,round(window_samples * (1. - overlap),digits=0))
    minind = 1:window_step:Nrows - window_samples
    overlap_factor = T(round(1. - overlap,digits=3))

    # allocate out array
    Aout = zeros(T,size(A))

    # loop through each window
    for ii in eachindex(minind)
        Ain = taper(A[minind[ii]:minind[ii]+window_samples-1,:],fs,
                    max_percentage=T(overlap_factor/2))
        Aout[minind[ii]:minind[ii]+window_samples-1,:] .+= ACF_kernel(@view(Ain[:,:]),g=g) .* overlap_factor
    end

    # windows at right edge
    reverseind = Nrows:-window_step:Nrows-window_samples
    for ii in eachindex(reverseind)
        Ain = taper(A[reverseind[ii]-window_samples+1:reverseind[ii],:],fs,
                    max_percentage=T(overlap_factor/2))
        Aout[reverseind[ii]-window_samples+1:reverseind[ii],:] .+= ACF_kernel(@view(Ain[:,:]),g=g) .* overlap_factor
    end

    copyto!(A,Aout)
    return nothing
end

adaptive_filter(A::AbstractArray{T}, window::AbstractFloat, fs::Real; g::AbstractFloat=2.,
                overlap::AbstractFloat=0.9) where T <: AbstractFloat = (U = deepcopy(A);
                adaptive_filter!(U,window,fs,g=g,overlap=overlap);
                return U)

adaptive_filter!(C::CorrData, window::AbstractFloat; g::AbstractFloat=2.,
                 overlap::AbstractFloat=0.9) = (adaptive_filter!(C.corr,
                 window,C.fs,g=g,overlap=overlap); return nothing)
adaptive_filter(C::CorrData, window::AbstractFloat; g::AbstractFloat=2.,
                overlap::AbstractFloat=0.9) = (U = deepcopy(C);
                adaptive_filter!(U.corr,window,U.fs,g=g,overlap=overlap);
                return U)

function ACF_kernel(A::AbstractArray{T}; g::AbstractFloat=2.) where T <: AbstractFloat
    Nrows, Ncols = size(A)
    g = T(g)
    # fft the 2D array
    spec = rfft(A,1)

    # create auto-covariance function
    Nspec = size(spec,1)
    S1 = zeros(complex(T),Nspec)
    S2 = zeros(complex(T),Nspec)
    for ii = 1:Ncols
        for jj = 1:Nspec
            S2[jj] += spec[jj,ii] .* conj(spec[jj,ii])
        end
    end

    for ii = 1:Ncols
        for jj = 1:Ncols
            if jj != ii
                for kk = 1:Nspec
                    S1[kk] += spec[kk,ii] .* conj(spec[kk,jj])
                end
            end
        end
    end

    # construct filter p
    p = ((S1 .- S2) ./ (S2 .* (Ncols -1))).^g

    # make ifft
    spec .*= p
    Aout = irfft(spec,Nrows,1)
    return Aout
end

"""

  robustpws(A)


Performs combined robust and phase-weighted stack on time series `A`.

Follows methods of Pavlis and Vernon, 2010 and Schimmel and Paulssen, 1997.
This method improves stacking by using the global convergence of the robust
stack and the phase convergence of the phase-weighted stack.
The robust stack downweights outlier correlations, while the phase-weighted
stack downweights outlier phases.

# Arguments
- `A::AbstractArray`: Time series stored in columns.
- `ϵ::AbstractFloat`: Threshold for convergence of robust stack.
- `maxiter::Int`: Maximum number of iterations to converge to robust stack.
- `pow::Int`: Sharpness of transition from phase similarity to dissimilarity.
"""
function robustpws(A::AbstractArray{T};ϵ::AbstractFloat=Float32(1e-6),
                     maxiter::Int=10,pow::Real=2) where T <: AbstractFloat
    N = size(A,2)
    Bold = median(A,dims=2)
    w = Array{T}(undef,N)
    r = Array{T}(undef,N)
    d2 = Array{T}(undef,N)

    # do 2-norm for all columns in A
    for ii = 1:N
        d2[ii] = norm(A[:,ii],2)
    end

    BdotD = sum(A .* Bold,dims=1)

    for ii = 1:N
        r[ii] = norm(A[:,ii] .- (BdotD[ii] .* Bold),2)
        w[ii] = abs(BdotD[ii]) ./ d2[ii] ./ r[ii]
    end
    w ./= sum(w)

    Bnew = mean(A,weights(w),dims=2)

    # check convergence
    ϵN = norm(Bnew .- Bold,2) / (norm(Bnew,2) * N)
    Bold = Bnew
    iter = 0
    while (ϵN > ϵ) && (iter <= maxiter)
        BdotD = sum(A .* Bold,dims=1)

        for ii = 1:N
            r[ii] = norm(A[:,ii] .- (BdotD[ii] .* Bold),2)
            w[ii] = abs(BdotD[ii]) ./ d2[ii] ./ r[ii]
        end
        w ./= sum(w)

        Bnew = mean(A,weights(w),dims=2)

        # check convergence
        ϵN = norm(Bnew .- Bold,2) / (norm(Bnew,2) * N)
        Bold = Bnew
        iter += 1
    end

    # multiply by weights
    W = A .* w'

    # return the phase weighted stack
    return pws(W,pow=pow)
end

robustpws!(C::CorrData; ϵ::AbstractFloat=Float32(1e-6),
           maxiter::Int=10,pow::Int=2) = (C.corr=robustpws(C.corr,ϵ=ϵ,
           maxiter=maxiter,pow=pow); C.t = C.t[1:1]; return nothing)
robustpws(C::CorrData; ϵ::AbstractFloat=Float32(1e-6),
           maxiter::Int=10,pow::Int=2) = (U = deepcopy(C);
           U.corr=robustpws(U.corr,ϵ=ϵ,
           maxiter=maxiter,pow=pow); U.t = U.t[1:1];return U)
