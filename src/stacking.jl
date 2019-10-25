export stack, stack!, remove_nan, remove_nan!, shorten, shorten!
export pws, pws!, robuststack, robuststack!
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
- `phase_smoothing::Float64`: Enables phase-weighted stacking. `phase_smoothing`
                              is the time window in seconds for phase smoothing
                              in the phase-weighted stack.
"""
function stack!(C::CorrData; interval::Union{Month,Day,Hour,Second}=Day(1),
                allstack::Bool=false,phase_smoothing::Float64=0.)

    if allstack == true
          if phase_smoothing == 0
                C.corr = mean(C.corr,dims=2)
          else # use phase-weighted stack
                C.corr = pws(C.corr, C.fs, power=2, timegate=phase_smoothing)
          end
          C.t = [C.t[1]]
    else # stack by interval
          C.t = d2u.(round.(u2d.(C.t),interval,RoundDown))
          stackT = unique(C.t)
          ind = indexin(stackT,C.t)
          push!(ind,length(C.t)+1)
          stack_out = Array{eltype(C.corr)}(undef,size(C.corr)[1],length(stackT))

          for ii = 1:length(stackT)
                 if phase_smoothing == 0
                       stack_out[:,ii] = mean(C.corr[:,ind[ii]:ind[ii+1]-1],dims=2)
                 else
                       stack_out[:,ii] = pws(C.corr[:,ind[ii]:ind[ii+1]-1],
                                             C.fs, power=2,
                                             timegate=phase_smoothing)
                 end
          end

          C.t = stackT
          C.corr = stack_out
    end
    return C
end
stack(C::CorrData; interval::Union{Month,Day,Hour,Second}=Day(1),
      allstack::Bool=false,phase_smoothing::Float64=0.) = (U = deepcopy(C);
      stack!(U,interval=interval,allstack=allstack,
             phase_smoothing=phase_smoothing);return U)

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

  pws(A,fs,power,timegate)

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
- `fs::Float64`: Sampling rate of time series (in Hz).
- `power::Int`: Sharpness of transition from phase similarity to dissimilarity.
- `timegate::Float64`: Phase stack smoothing window (in seconds).
"""
function pws(A::AbstractArray, fs::Float64; power::Int=2, timegate::Float64=1.)
      M,N = size(A)
      analytic = A .+ im .* hilbert(A)
      phase = angle.(analytic)
      phase_stack = mean(exp.(im.*phase),dims=2)[:,1] ./ N # reduce array dimension
      phase_stack = abs.(phase_stack).^power

      # smoothing
      timegate_samples = convert(Int,timegate * fs * 0.5)
      phase_stack = movingaverage(phase_stack,timegate_samples)
      abs_max!(phase_stack)
      weighted = mean(A .* phase_stack,dims=2)
      return weighted
end
pws!(C::CorrData; power::Int=2, timegate::AbstractFloat=1.) = (C.corr=pws(C.corr,
     C.fs,power=power,timegate=timegate); C.t = C.t[1:1]; return C)
pws(C::CorrData; power::Int=2, timegate::AbstractFloat=1.) = (U = deepcopy(C);
    U.corr=pws(U.corr,U.fs,power=power,timegate=timegate); U.t = U.t[1:1];
    return U)

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
