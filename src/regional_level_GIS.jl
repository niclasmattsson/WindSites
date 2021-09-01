function plotdist(gisregion::String, minturbine=1000, maxpopdens=150; args...)
    filename = "regionalwindGIS_$(gisregion)_minturbine=$(minturbine)_maxpopdens=$(maxpopdens).csv"
    df = DataFrame(CSV.File(in_datafolder("output", filename)))
    plotdist(df; title=gisregion, args...)
end

function plotdist(gisregions::Vector{String}, minturbine=1000, maxpopdens=150; args...)
    dfs = [DataFrame(CSV.File(in_datafolder("output", "regionalwindGIS_$(gisregion)_minturbine=$(minturbine)_maxpopdens=$(maxpopdens).csv")))
            for gisregion in gisregions]
    for (i, gisregion) in enumerate(gisregions)
        dfs[i][:,:gisregion] .= i
    end
    df = vcat(dfs...)
    plotdist(df; title=join(gisregions, ", "), args...)
end

function plotdist(df::DataFrame; mincapac=0, area_mult=1.0, scatterplot=0, bins=100,
                    variable=:exploit_tot, alpha=1, markersize=2, line=:mean,
                    scale=1, xmax=30, title="", legend=true)
    df = df[df[!,:capac].>=mincapac, :]
    plotly()
    blank = "<span style='color:white; font-size: 1px;'>q</span>"
    title = replace(title, "GADM" => "")
    title = replace(title, r"\d" => "")
    if scatterplot == 0
        p = histogram(df[!,variable], bins=range(0, xmax+0.5, length=bins),
            xlims=(0,Inf), ylims=(0,Inf), title=title,
            xlabel="Exploited area per municipality/county [%]",
            ylabel="Number of municipalities/counties<br>$blank",
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            size=scale.*(800,550), legend=false)  # title=variable
        line_x = line == :mean ? mean(df[!,variable]) : 100*area_mult
        line != :none && plot!(line_x*[1, 1], collect(ylims(p)), line=(3, :dash))
        display(p)
        return df
    end
    # Filter out NaNs so scatter plots don't get messed up.
    # https://github.com/JuliaPlots/Plots.jl/issues/3258
    df = df[.!(isnan.(df[!,:capac]) .| isnan.(df[!,:class]) .| isnan.(df[!,variable])), :]
    if scatterplot == 1
        p = scatter(df[!,variable], df[!,:capac], 
            xlabel="Exploited area per municipality/county [%]",
            ylabel="MW", title=title,
            size=scale.*(800,550), legend=false, colorbar=true,
            markersize=scale*df[!,:class].^(1+markersize/10),
            marker_z=df[!,:masked], color=:watermelon, alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable],digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        display(p)
    elseif scatterplot == 2
        p = scatter(df[!,variable], df[!,:class], 
            xlabel="Exploited area per municipality/county [%]", xlims=(0, xmax+0.5),
            ylabel="Mean wind class",
            size=scale.*(800,550), legend=false, colorbar=true,
            markersize=scale*df[!,:capac].^(markersize/7),
            markerstrokecolor=RGBA(0,0,0,.2), markerstrokewidth=1,
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            marker_z=df[!,:masked], color=:plasma, alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable],digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        display(p)
    elseif scatterplot == 3
        p = scatter(df[!,variable], df[!,:class], 
            xlabel="Exploited area per municipality/county [%]", xlims=(0, xmax+0.5),
            ylabel="Mean wind class", title=title,
            size=scale.*(800,550), legend=false, colorbar=true,
            markersize=markersize*scale,
            markerstrokecolor=RGBA(0,0,0,.2), markerstrokewidth=0,
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            color=RGB(0,0,0), alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable],digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        display(p)
    elseif scatterplot == 4
        cc = [RGB(0,0,1), RGB(0,0,0), RGB(0,1,0), RGB(1,0,0)]
        p = scatter(df[!,variable], df[!,:windspeed], 
            xlabel="Exploited area per municipality/county [%]",
            ylabel="Mean wind speed [m/s]", xlims=(0, xmax+0.5), ylims=(3.5, 11), title=title,
            size=scale.*(950,550), label="", colorbar=true,
            markersize=markersize*scale,
            markerstrokecolor=RGBA(0,0,0,.2), markerstrokewidth=0,
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            color=cc[df[!,:gisregion]], alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable],digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        qq = [-1 -1 -1 -1]
        scatter!(qq,qq,color=cc[[1 2 3 4]],alpha=alpha+0.1,legend=:outertopright,
            markerstrokewidth=0, legendfont=12*scale, label=permutedims(split(title, ", ")))    
        display(p)
    elseif scatterplot == 5
        cc = [RGB(0,0,1), RGB(0,0,0), RGB(0,1,0), RGB(1,0,0)]
        p = scatter(df[!,variable], df[!,:windspeed], 
            xlabel="Exploited area per municipality/county [%]",
            ylabel="Mean wind speed [m/s]", xlims=(0, xmax+0.5), ylims=(3.5, 11), title=title,
            size=scale.*(800+150*legend,550), label="", colorbar=false,
            markersize=scale*df[!,:capac].^(markersize/7),
            markerstrokecolor=RGBA(0,0,0,.2), markerstrokewidth=0,
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            color=cc[df[!,:gisregion]], alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable],digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        qq = [-1 -1 -1 -1]
        legend && scatter!(qq,qq,color=cc[[1 2 3 4]],alpha=alpha+0.1,legend=:outertopright,
            markerstrokewidth=0, legendfont=12*scale, label=permutedims(split(title, ", ")))
        display(p)
    elseif scatterplot == 6
        cc = [RGB(0,0,1), RGB(0,0,0), RGB(0,1,0), RGB(1,0,0)]
        p = scatter(df[!,variable], max.(-1, log10.(df[!,:popdens])), 
            xlabel="Exploited area per municipality/county [%]",
            ylabel="Population density [log10 persons/km2]", xlims=(0, xmax+0.5), ylims=(-1.1, Inf), title=title,
            size=scale.*(950,550), label="", colorbar=false,
            markersize=scale*df[!,:capac].^(markersize/7),
            markerstrokecolor=RGBA(0,0,0,.2), markerstrokewidth=0,
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            color=cc[df[!,:gisregion]], alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable],digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        qq = [-1 -1 -1 -1]
        scatter!(qq,qq,color=cc[[1 2 3 4]],alpha=alpha+0.1,legend=:outertopright,
            markerstrokewidth=0, legendfont=12*scale, label=permutedims(split(title, ", ")))
        display(p)
    end 
    return df
end

# function scatterplots_model(gisregion, type=[:total]; showlines=false)
#     dfg = DataFrame(CSV.File(in_datafolder("output", "regionalwindGIS_$gisregion.csv")))
#     dfm = DataFrame(CSV.File(in_datafolder("output", "windresults_$gisregion.csv")))
#     dfm[!,:mcap] = vec(sum(Array(dfm[:, [:cap1, :cap2, :cap3, :cap4, :cap5]]), dims=2))
#     dfm[!,:moffcap] = vec(sum(Array(dfm[:, [:offcap1, :offcap2, :offcap3, :offcap4, :offcap5]]), dims=2))
#     dfm[!,:melec] = vec(sum(Array(dfm[:, [:elec1, :elec2, :elec3, :elec4, :elec5]]), dims=2))
#     dfm[!,:moffelec] = vec(sum(Array(dfm[:, [:offelec1, :offelec2, :offelec3, :offelec4, :offelec5]]), dims=2))
#     dfm[:mclass] = (1*dfm[:cap1] + 2*dfm[:cap2] + 3*dfm[:cap3] + 4*dfm[:cap4] + 5*dfm[:cap5]) ./ dfm[:mcap]
#     dfm[:mclass] = replace(round.(dfm[:mclass], digits=2), NaN => missing)
#     dfm[:moffclass] = (1*dfm[:offcap1] + 2*dfm[:offcap2] + 3*dfm[:offcap3] + 4*dfm[:offcap4] + 5*dfm[:offcap5]) ./ dfm[:moffcap]
#     dfm[:moffclass] = replace(round.(dfm[:moffclass], digits=2), NaN => missing)

#     distancevars = matread(in_datafolder("output", "distances_$gisregion.mat"))
#     connectedregions = findall(distancevars["connected"] .| distancevars["connectedoffshore"])
#     connected = [c for c in connectedregions if c[1] < c[2]]

#     df = innerjoin(dfg, dfm, on=:region)
#     select!(df, [:region, :capac, :offcapac, :elec2018, :class, :mcap, :moffcap,
#                     :melec, :moffelec, :mclass, :moffclass, :exploit_tot, :masked])
#     plotly()
#     println(type)
#     if :onshore in type
#         s = scatter(df[:capac], df[:mcap], markersize=df[:mclass]*2, title="Onshore capacity",
#             hover=df[:region].*" class ".*string.(df[:mclass]), legend=false, colorbar=true,
#             xlabel="GW (real)", ylabel="GW (model)", size=(800,550))
#         display(s)
#     elseif :offshore in type
#         s = scatter(df[:offcapac], df[:moffcap], markersize=df[:moffclass]*2, title="Offshore capacity",
#             hover=df[:region].*" class ".*string.(df[:moffclass]), legend=false, colorbar=true,
#             xlabel="GW (real)", ylabel="GW (model)", size=(800,550))
#         display(s)
#     elseif :total in type
#         x = df[:capac] + df[:offcapac]
#         y = df[:mcap] + df[:moffcap]
#         s = plot()
#         xx = [x[c[d]] for d in 1:2, c in connected]
#         yy = [y[c[d]] for d in 1:2, c in connected]
#         showlines && plot!(xx, yy, c=WindSites.RGB(.85,.85,.85), legend=:none)
#         scatter!(x, y, markersize=df[:class]*2, title="Total capacity", c=1,
#             hover=df[:region].*"<br>real onshore = ".*string.(df[:capac]).*" MW".*
#             "<br>real offshore = ".*string.(df[:offcapac]).*" MW".*
#             "<br>exploited onshore = ".*string.(df[:exploit_tot]).*"%".*
#             "<br>class = ".*string.(df[:class]).*
#             "<br>m_onshore = ".*string.(df[:mcap]).*" GW<br>m_offshore = ".*string.(df[:moffcap]).*
#             " GW<br>m_class = ".*string.(df[:mclass]).*
#             "<br>m_offclass = ".*string.(df[:moffclass]).*"<br>masked onshore = ".*string.(df[:masked]),
#             legend=false, colorbar=true, marker_z=df[:masked], color=:watermelon, 
#             xlabel="GW (real)", ylabel="GW (model)", size=(800,550))
#         display(s)
#     elseif :elec in type
#         s = scatter(df[:elec2018], df[:melec], markersize=df[:mclass], xlabel="GWh (real)", ylabel="GWh (model)")
#         display(s)
#     end
#     df
#  end
