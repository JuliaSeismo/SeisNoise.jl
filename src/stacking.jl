export stack, stack!
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
