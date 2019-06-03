export plot

"""

  plot(C)

Plot a CorrData matrix.
"""
function plot(C::CorrData;kwargs...)
    lags = -C.maxlag:1/C.fs:C.maxlag
    times = timestamp.(C.t)
    heatmap(lags,times,C.corr',c=:RdBu,kwargs...)
    xlabel!("Lags [s]")
end
