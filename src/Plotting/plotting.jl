@recipe function f(C::CorrData)
    # set up the subplots
    grid := false
    layout := (2,1)

    lags = -C.maxlag:1/C.fs:C.maxlag
    times = Dates.format.(Dates.unix2datetime.(C.t),"yyyy/m/d HH:MM")

    # main heatmap
    @series begin
        legend := false
        seriestype := :heatmap
        ytickfontsize --> 10 
        xtickfontsize --> 8
        subplot := 1
        seriescolor --> :balance
        lags,times,C.corr'
    end

    # bottom stack
    @series begin
        subplot := 2
        label := false
        title --> C.name
        xguide --> "Lag [s]"
        titlefontsize --> 6 
        ytickfontsize --> 10 
        xtickfontsize --> 8
        seriescolor --> :black
        xlims --> (lags[1],lags[end])
        lags,stack(C,allstack=true).corr
    end
end