export corrplot

"""

  corrplot(C)

Plot a CorrData matrix.
"""
function corrplot(C::CorrData)
    lags = -C.maxlag:1/C.fs:C.maxlag
    times = Dates.format.(Dates.unix2datetime.(C.t),"yyyy/m/d HH:MM")
    Cstack = stack(C,allstack=true)
    plot(
         heatmap(lags,times,C.corr',c=:balance,legend=:none),
         plot(lags,Cstack.corr,c=:black,linewidth=1.5,legend=:none,xlabel="Lag [s]"),
         layout = grid(2,1,heights=[0.75,0.25]), link=:x,dpi=1000)
end
