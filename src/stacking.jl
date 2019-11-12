export stack, stack!, remove_nan, remove_nan!, shorten, shorten!
export pws, pws!, robuststack, robuststack!, adaptive_filter!, adaptive_filter
"""

  stack!(C)

Stack correlation by time interval.
The default is to stack by day. Using `allstack == true` will stack all available
correlations. To use phase-weighted stack, specify the amount of
`phase_smoothing` in seconds.

# Arguments
- `C::CorrData`: Correlation data.
- `interval::Union{Month,Day,Hour,Second}`: Interval over which to stack `C`.
- `allstack::Bool`: If `true`, stack all data.
- `stack_type::Function`: Type of stacking. Options are mean, pws, robuststack, etc..
"""
function stack!(C::CorrData; interval::Union{Month,Day,Hour,Second}=Day(1),
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
          stack_out = Array{eltype(C.corr)}(undef,size(C.corr)[1],length(stackT))

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
stack(C::CorrData; interval::Union{Month,Day,Hour,Second}=Day(1),
      allstack::Bool=false,stacktype::Function=mean) = (U = deepcopy(C);
      stack!(U,interval=interval,allstack=allstack,
             stacktype=stacktype);return U)

"""

  robuststack(A,ϵ)

Performs robust stack on array `A`.

Follows methods of Pavlis and Vernon, 2010.

# Arguments
- `A::AbstractArray`: Time series stored in columns.
- `ϵ::AbstractFloat`: Threshold for convergence of robust stack.
"""
function robuststack(A::AbstractArray{T};ϵ::AbstractFloat=eps(Float32)) where T <: AbstractFloat
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
    while ϵN < ϵ
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
    return Bnew
end
robuststack!(C::CorrData;ϵ::AbstractFloat=eps(Float32)) =
       (C.corr = robuststack(C.corr,ϵ=ϵ); C.t = C.t[1:1]; return C)
robuststack(C::CorrData;ϵ::AbstractFloat=eps(Float32))  =
       (U = deepcopy(C); U.corr = robuststack(U.corr,ϵ=ϵ); U.t = U.t[1:1];
       return U)

"""

  pws(A,power=2)

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
- `power::Int`: Sharpness of transition from phase similarity to dissimilarity.
"""
function pws(A::AbstractArray; power::Int=2)
      M,N = size(A)
      phase_stack = abs.(mean(exp.(im.*angle.(A .+ im .* hilbert(A))),dims=2)[:,1] ./ N).^power
      abs_max!(phase_stack)
      return mean(A .* phase_stack,dims=2)
end
pws!(C::CorrData; power::Int=2) = (C.corr=pws(C.corr,power=power); C.t = C.t[1:1]; return nothing)
pws(C::CorrData; power::Int=2) = (U = deepcopy(C);
    U.corr=pws(U.corr,power=power); U.t = U.t[1:1];return U)

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
      C.corr = C.corr[:,ind]
      C.t = C.t[ind]
      return nothing
end
remove_nan(C::CorrData) = (U = deepcopy(C);remove_nan!(U);return U)

"""

  shorten!(C)

Clip CorrData `C` from lags τ = 0 to abs(maxlag).
"""
function shorten!(C::CorrData, maxlag::Float64)
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
shorten(C::CorrData,maxlag::Float64) = (U = deepcopy(C); shorten!(U,maxlag);return U)

"""

  adaptive_filter!(A)

Adaptive covariance filter to enhance coherent signals.

Fellows the method of Nakata et al., 2015 (Appendix B). The filtered signal
`X` is given by `X = irfft(P*X(w))`` where `X` is the FFT'ed spectra of `A`
and `P` is the filter. `P` is constructed by using the temporal covariance matrix.
​
# Arguments
- `A::AbstractArray`: Array of daily/hourly cross-correlation functions
- `g::Int`: Positive number to adjust the filter harshness
"""
function adaptive_filter!(A::AbstractArray{T}; g::Int=2) where T <: AbstractFloat
    if ndims(A) == 1
        return A
    end

    Nrows, Ncols = size(A)

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
    A .= irfft(spec,Nrows,1)
    return nothing
end
adaptive_filter(A::AbstractArray,g::Int) = (U = deepcopy(A);adaptive_filter!(U,g=g);
                return U)
