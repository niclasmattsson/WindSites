module WindSites

using Proj4, XLSX, CSV, DataFrames, Dates, GDAL_jll, Pkg.TOML, Plots, MAT,
        Loess, SmoothingSplines, Plots.PlotMeasures

export openmap, readusa, readdk, readuk, readse, readde, readturbinedata, shapefile2csv,
    scatterplots_model, plotdist, fill_missing_rotordiams!, closestturbine, turbines2parks

function openmap(df::DataFrame, turbinenumber::Int)
    openmap(df, turbinenumber, :bing)
    openmap(df, turbinenumber, :google)
end
    
function openmap(df::DataFrame, turbinenumber::Int, source::Symbol)
    lon, lat = df[turbinenumber, :lon], df[turbinenumber, :lat]
    if source == :google
        url = "http://maps.google.com/maps?t=k\"&\"q=loc:$lat+$lon"
    else
        url = "https://bing.com/maps/default.aspx?cp=$lat~$lon\"&\"lvl=18\"&\"style=a\"&\"sp=point.$(lat)_$(lon)_"
    end
    c = Cmd(`cmd /c start \"\" $url`, windows_verbatim=true)
    run(c)
    return df[turbinenumber, :]
end

function openmap(lon::Real, lat::Real)
    url = "http://maps.google.com/maps?t=k\"&\"q=loc:$lat+$lon"
    # Extra quotes to avoid errors with special chars ? and &:
    # https://superuser.com/questions/36728/can-i-launch-urls-from-command-line-in-windows
    # https://discourse.julialang.org/t/quoting-special-characters-of-a-url-in-cmd-objects-on-windows/44324
    c = Cmd(`cmd /c start \"\" $url`, windows_verbatim=true)
    run(c)
end

function openmapchrome(lon::Real, lat::Real)
    command = "C:/Program Files (x86)/Google/Chrome/Application/chrome"
    url = "http://maps.google.com/maps?t=k&q=loc:$lat+$lon"
    println("Opening $url...")
    run(`$command $url`)
end

function readusa()
    # US Wind Turbine Database
    # https://eerscmap.usgs.gov/uswtdb/
    datadir = joinpath(@__DIR__, "..", "data")
    df = DataFrame!(CSV.File("$datadir/USA_uswtdb_v3_1_20200717.csv"))
    select!(df, [:t_cap, :xlong, :ylat, :p_year, :t_model, :t_hh, :t_rd])
    rename!(df, [:capac, :lon, :lat, :year, :model, :hubheight, :rotordiam])
    delete!(df, ismissing.(df[:, :capac]))
    delete!(df, ismissing.(df[:, :year]))
    df[!, :capac] = convert.(Int, df[!, :capac])
    return df
end

parse_not_missing(T, x) = ismissing(x) ? missing : parse(T, x)
convert_not_missing(T, x) = ismissing(x) ? missing : convert(T, x)

function readuk()
    # UK Renewable Energy Planning Database (REPD)
    # https://www.gov.uk/government/publications/renewable-energy-planning-database-monthly-extract
    datadir = joinpath(@__DIR__, "..", "data")
    cols = [6,9,14,15,16,19,25,26,47]
    df = DataFrame!(CSV.File("$datadir/UK_REPD-June-2020-update.csv"))[:, cols]
    rename!(df, [:type, :capac_park, :capac, :nturbines, :height, :status, :lon, :lat, :year])
    delete!(df, .!startswith.(df[!, :type], "Wind"))
    delete!(df, df[:, :status] .!= "Operational")
    select!(df, Not(:status))
    df[!, :type] = replace.(df[!, :type], "Wind " => "")
    df[!, :nturbines] = parse.(Int, df[!, :nturbines])
    df[!, :height] = parse_not_missing.(Int, df[!, :height])
    df[!, :capac] = parse.(Float64, df[!, :capac])
    df[!, :capac_park] = parse.(Float64, replace.(df[!, :capac_park], "," => ""))
    df[!, :lon] = parse.(Int, replace.(df[!, :lon], "," => ""))
    df[!, :lat] = parse.(Int, replace.(df[!, :lat], "," => ""))
    df[!, :year] = Date.(df[!, :year], dateformat"d/m/y")

    x = df[!, :lon] # BNG: British National Grid = EPSG:27700
    y = df[!, :lat] # https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid#General
    source = Projection("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs")
    dest = Projection("+proj=longlat +datum=WGS84 +no_defs")
    len = size(df,1)
    lon, lat = zeros(len), zeros(len)
    for i = 1:len
        lon[i], lat[i], _ = Proj4.transform(source, dest, [x[i] y[i] 0])
    end
    df[!, :lon] = lon
    df[!, :lat] = lat
    return df
end

function readdk()
    # DK Energistyrelsen: Stamdataregister for vindkraftanlæg
    # https://ens.dk/service/statistik-data-noegletal-og-kort/data-oversigt-over-energisektoren
    datadir = joinpath(@__DIR__, "..", "data")
    cols = [3,13,14,2,7,5,4,10,58]
    df = DataFrame!(XLSX.readtable("$datadir/DK_anlaegprodtilnettet_0.xlsx", "IkkeAfmeldte-Existing turbines",
        "A:BM", first_row=19, header=false, infer_eltypes=false)...)[!, cols]
    rename!(df, [:capac, :lon, :lat, :year, :model, :hubheight, :rotordiam, :onshore, :elec2018])
    df[!, :capac] = convert.(Int, round.(df[!, :capac]))
    df[!, :model] = string.(df[:, :model])
    df[!, :hubheight] = convert.(Float64, df[!, :hubheight])
    df[!, :rotordiam] = convert.(Float64, df[!, :rotordiam])
    df[!, :year] = Dates.year.(df[!, :year])
    df[!, :onshore] = uppercase.(df[!, :onshore]) .== "LAND"    # no missings here, no need to check
    df[!, :elec2018] = df[!, :elec2018] / 1000
    delete!(df, ismissing.(df[:, :lon]))

    x = df[!, :lon] # European Terrestrial Reference System 1989: ETRS89/UTM zone 32N = EPSG:25832
    y = df[!, :lat] # https://en.wikipedia.org/wiki/European_Terrestrial_Reference_System_1989
    source = Projection("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    dest = Projection("+proj=longlat +datum=WGS84 +no_defs")
    len = size(df,1)
    lon, lat = zeros(len), zeros(len)
    for i = 1:len
        lon[i], lat[i], _ = Proj4.transform(source, dest, [x[i] y[i] 0])
    end
    df[!, :lon] = lon
    df[!, :lat] = lat
    return df
end

function xlsx2csv()
    datadir = joinpath(@__DIR__, "..", "data")
    df = DataFrame!(XLSX.readtable("$datadir/SE_Vindbrukskollen_export_allman_Prod.xlsx", "Vindkraftverk",
        "A:AN", first_row=1, header=true, infer_eltypes=true)...)
    CSV.write("$datadir/SE_Vindbrukskollen_export_allman_Prod.csv", df)
    return df
end

function readse()
    # Länsstyrelsen: Vindbrukskollen
    # https://vbk.lansstyrelsen.se/  (click "Excel-export")
    # shapefile: https://ext-dokument.lansstyrelsen.se/gemensamt/geodata/ShapeExport/LST.VKOLLEN_VINDKRAFTVERK.zip 
    # shapefile2csv("LST.VKOLLEN_VINDKRAFTVERK/LST_VKOLLEN_VVERK.shp")
    # but actually used xlsx2csv() above to convert the Excel file 
    datadir = joinpath(@__DIR__, "..", "data")
    df = DataFrame(CSV.File("$datadir/SE_Vindbrukskollen_export_allman_Prod.csv"))
    select!(df, ["Status", "Placering", "E-Koordinat", "N-Koordinat", "Totalhöjd (m)", "Navhöjd (m)",
        "Rotordiameter (m)", "Maxeffekt (MW)", "Uppfört", "Fabrikat", "Modell"])
    rename!(df, [:status, :type, :lon, :lat, :height, :hubheight, :rotordiam, :capac, :year, :brand, :model])
    delete!(df, df[:, :status] .!= "Uppfört")
    select!(df, Not(:status))
    twoturbines = findall(ismissing.(df[!, :type]))
    if twoturbines == [1196, 1197]     # hard code a fix for two turbines
        df[twoturbines, :type] .= "Land"
    end
    df[!, :onshore] = startswith.(df[!, :type], "Land")     # all turbines are marked "Land" or "Vatten"
    select!(df, Not(:type))
    delete!(df, ismissing.(df[:, :capac]))
    delete!(df, df[:, :capac] .== 0)
    df[!, :capac] = Int.(df[!, :capac] * 1000)
    # df[!, :year] = [d == "Jan 01, 1900" ? missing : Date(d, dateformat"u d, y") for d in df[!, :year]]
    df[df[!, :year] .< Date(1980), :year] .= missing
    df[!, :model] = strip.(coalesce.(df[!, :brand], "")) .* " " .* strip.(coalesce.(df[!, :model], ""))
    df[!, :model] = [m == " " ? missing : m for m in df[!, :model]]
    select!(df, Not(:brand))

    x = df[!, :lon] # SWEREF99 TM = EPSG:3006
    y = df[!, :lat] # https://www.lantmateriet.se/sv/Kartor-och-geografisk-information/gps-geodesi-och-swepos/referenssystem/tvadimensionella-system/sweref-99-projektioner/
    source = Projection("+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    dest = Projection("+proj=longlat +datum=WGS84 +no_defs")
    len = size(df,1)
    lon, lat = zeros(len), zeros(len)
    for i = 1:len
        lon[i], lat[i], _ = Proj4.transform(source, dest, [x[i] y[i] 0])
    end
    df[!, :lon] = lon
    df[!, :lat] = lat
    return df
end

function readturbinedata(; rotordiamplots=false)
    df_dk = readdk()
    df_usa = readusa()
    df_uk = readuk()
    df_se = readse()
    df_de = readde()

    fill_missing_rotordiams!(df_dk, df_usa, df_uk, df_se, df_de; rotordiamplots=rotordiamplots)

    CSV.write(in_datafolder("turbines_DK.csv"), df_dk)
    CSV.write(in_datafolder("turbines_USA.csv"), df_usa)
    CSV.write(in_datafolder("turbines_UK.csv"), df_uk)
    CSV.write(in_datafolder("turbines_SE.csv"), df_se)
    CSV.write(in_datafolder("turbines_DE.csv"), df_de)

    return df_dk, df_usa, df_uk, df_se, df_de
end

function readde()
    # Eichhorn et al. (2019), Spatial Distribution of Wind Turbines, Photovoltaic Field Systems, Bioenergy,
    # and River Hydro Power Plants in Germany
    # https://www.mdpi.com/2306-5729/4/1/29
    # https://www.ufz.de/record/dmp/archive/5467/en/

    # shapefile2csv("renewable_energy_plants_germany_until_2015/windpower.shp")
    datadir = joinpath(@__DIR__, "..", "data")
    df = DataFrame!(CSV.File("$datadir/DE_windpower.csv"))[!, 2:7]
    rename!(df, [:year, :capac, :hubheight, :rotordiam, :lon, :lat])

    x = df[!, :lon] # European Terrestrial Reference System 1989: ETRS89/LAEA Europe = EPSG:3035
    y = df[!, :lat] # https://en.wikipedia.org/wiki/European_Terrestrial_Reference_System_1989
    source = Projection("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    dest = Projection("+proj=longlat +datum=WGS84 +no_defs")
    len = size(df,1)
    lon, lat = zeros(len), zeros(len)
    for i = 1:len
        lon[i], lat[i], _ = Proj4.transform(source, dest, [x[i] y[i] 0])
    end
    df[!, :lon] = lon
    df[!, :lat] = lat
    return df
end

function shapefile2csv(filename)
    datadir = joinpath(@__DIR__, "..", "data")
    shapefile = "$datadir/$filename"
    layername = splitext(basename(filename))[1]
    # sql = "select FID+1 as FID from \"$layername\""

    ogrinfo_path() do ogrinfo
        @time run(`$ogrinfo -al -so $shapefile`)
    end

    ogr2ogr_path() do ogr2ogr
        # @time run(`$ogr2ogr -f CSV $outfile -dialect SQlite -sql $sql $shapefile`)
        @time run(`$ogr2ogr -f CSV $datadir/$layername.csv $shapefile`)
    end
end

# println("\nProjecting coordinates (Mollweide)...")
# res = 0.01
# res2 = res/2
# lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
# lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
# source = LonLat()
# dest = Projection("+proj=moll +lon_0=$(mean(lons))")
# xs, ys = xygrid(lons, lats)
# Proj4.transform!(source, dest, vec(xs), vec(ys))

in_datafolder(names...) = joinpath(getconfig("datafolder"), names...)

getconfig(key) = getconfig()[key]

function getconfig()
    configfile = joinpath(homedir(), ".GlobalEnergyGIS_config")
    if !isfile(configfile)
        error("Configuration file missing, please run saveconfig(datafolder, uid, api_key) first. See GlobalEnergyGIS README.")
    end
    return TOML.parsefile(configfile)
end

# estimate missing rotor diameters from turbine capacity using smoothing splines
# fit separate splines to turbines below and above 1 MW
function fill_missing_rotordiams!(df_dk, df_usa, df_uk, df_se, df_de; rotordiamplots=false)
    df = df_de[:, [:rotordiam, :capac]]
    badrows = (df[!, :capac] .< 80) .& (df[!, :rotordiam] .> 20)
    df = df[.!badrows, :]
    for dd in [df_se, df_dk, df_usa]
        df = vcat(df, dd[:, [:rotordiam, :capac]])
    end
    rotordiamplots && plotly()
    rotordiamplots && display(histogram(df[ismissing.(df[!, :rotordiam]), :capac]))

    dfd = dropmissing(df, [:rotordiam, :capac])
    dfd = dfd[dfd[!, :rotordiam] .> 5, :]
    x = float.(dfd[:, :capac])
    y = float.(dfd[:, :rotordiam])
    rotordiamplots && scatter(x, y .+ 1*randn(size(dfd,1)), markersize=1, alpha=0.1)
    ss1 = fit(SmoothingSpline, x, y, 0.0001*maximum(x)^3)
    ss2 = fit(SmoothingSpline, x, y, 0.1*maximum(x)^3)
    # ss = loess(x, y, span=0.75, degree=3)
    xx = range(minimum(x), maximum(x), length=1000)
    yy1 = SmoothingSplines.predict(ss1, xx)
    yy2 = SmoothingSplines.predict(ss2, xx)
    rotordiamplots && display(plot!(xx, [yy1, yy2]))
    
    df_uk[:, :rotordiam] = zeros(Int, size(df_uk,1))
    for dd in [df_se, df_dk, df_de, df_uk, df_usa]
        rows = (dd[!, :capac] .<= 1000)
        is_uk = (size(dd,1) == size(df_uk,1)) 
        rows = is_uk ? rows : rows .& ismissing.(dd[!, :rotordiam])
        dd[rows, :rotordiam] =
            round.(Int, SmoothingSplines.predict(ss1, float.(dd[rows, :capac])))
        rows = (dd[!, :capac] .> 1000)
        rows = is_uk ? rows : rows .& ismissing.(dd[!, :rotordiam])
        dd[rows, :rotordiam] =
            round.(Int, SmoothingSplines.predict(ss2, float.(dd[rows, :capac])))
        dd[!, :rotordiam] = round.(Int, dd[:, :rotordiam])
    end
end

function scatterplots_model(gisregion, type=[:total]; showlines=false)
    dfg = DataFrame!(CSV.File(in_datafolder("output", "regionalwindGIS_$gisregion.csv")))
    dfm = DataFrame!(CSV.File(in_datafolder("output", "windresults_$gisregion.csv")))
    dfm[!,:mcap] = vec(sum(Array(dfm[:, [:cap1, :cap2, :cap3, :cap4, :cap5]]), dims=2))
    dfm[!,:moffcap] = vec(sum(Array(dfm[:, [:offcap1, :offcap2, :offcap3, :offcap4, :offcap5]]), dims=2))
    dfm[!,:melec] = vec(sum(Array(dfm[:, [:elec1, :elec2, :elec3, :elec4, :elec5]]), dims=2))
    dfm[!,:moffelec] = vec(sum(Array(dfm[:, [:offelec1, :offelec2, :offelec3, :offelec4, :offelec5]]), dims=2))
    dfm[:mclass] = (1*dfm[:cap1] + 2*dfm[:cap2] + 3*dfm[:cap3] + 4*dfm[:cap4] + 5*dfm[:cap5]) ./ dfm[:mcap]
    dfm[:mclass] = replace(round.(dfm[:mclass], digits=2), NaN => missing)
    dfm[:moffclass] = (1*dfm[:offcap1] + 2*dfm[:offcap2] + 3*dfm[:offcap3] + 4*dfm[:offcap4] + 5*dfm[:offcap5]) ./ dfm[:moffcap]
    dfm[:moffclass] = replace(round.(dfm[:moffclass], digits=2), NaN => missing)

    distancevars = matread(in_datafolder("output", "distances_$gisregion.mat"))
    connectedregions = findall(distancevars["connected"] .| distancevars["connectedoffshore"])
    connected = [c for c in connectedregions if c[1] < c[2]]

    df = innerjoin(dfg, dfm, on=:region)
    select!(df, [:region, :capac, :offcapac, :elec2018, :class, :mcap, :moffcap,
                    :melec, :moffelec, :mclass, :moffclass, :exploit_tot, :masked])
    plotly()
    println(type)
    if :onshore in type
        s = scatter(df[:capac], df[:mcap], markersize=df[:mclass]*2, title="Onshore capacity",
            hover=df[:region].*" class ".*string.(df[:mclass]), legend=false, colorbar=true,
            xlabel="GW (real)", ylabel="GW (model)", size=(800,550))
        display(s)
    elseif :offshore in type
        s = scatter(df[:offcapac], df[:moffcap], markersize=df[:moffclass]*2, title="Offshore capacity",
            hover=df[:region].*" class ".*string.(df[:moffclass]), legend=false, colorbar=true,
            xlabel="GW (real)", ylabel="GW (model)", size=(800,550))
        display(s)
    elseif :total in type
        x = df[:capac] + df[:offcapac]
        y = df[:mcap] + df[:moffcap]
        s = plot()
        xx = [x[c[d]] for d in 1:2, c in connected]
        yy = [y[c[d]] for d in 1:2, c in connected]
        showlines && plot!(xx, yy, c=WindSites.RGB(.85,.85,.85), legend=:none)
        scatter!(x, y, markersize=df[:class]*2, title="Total capacity", c=1,
            hover=df[:region].*"<br>real onshore = ".*string.(df[:capac]).*" MW".*
            "<br>real offshore = ".*string.(df[:offcapac]).*" MW".*
            "<br>exploited onshore = ".*string.(df[:exploit_tot]).*"%".*
            "<br>class = ".*string.(df[:class]).*
            "<br>m_onshore = ".*string.(df[:mcap]).*" GW<br>m_offshore = ".*string.(df[:moffcap]).*
            " GW<br>m_class = ".*string.(df[:mclass]).*
            "<br>m_offclass = ".*string.(df[:moffclass]).*"<br>masked onshore = ".*string.(df[:masked]),
            legend=false, colorbar=true, marker_z=df[:masked], color=:watermelon, 
            xlabel="GW (real)", ylabel="GW (model)", size=(800,550))
        display(s)
    elseif :elec in type
        s = scatter(df[:elec2018], df[:melec], markersize=df[:mclass], xlabel="GWh (real)", ylabel="GWh (model)")
        display(s)
    end
    df
 end

function plotdist(gisregion::String; args...)
    df = DataFrame!(CSV.File(in_datafolder("output", "regionalwindGIS_$gisregion.csv")))
    plotdist(df; args...)
end

function plotdist(gisregions::Vector{String}; args...)
    dfs = [DataFrame!(CSV.File(in_datafolder("output", "regionalwindGIS_$gisregion.csv")))
            for gisregion in gisregions]
    df = vcat(dfs...)
    plotdist(df; args...)
end

function plotdist(df::DataFrame; mincapac=0, area_mult=1.0, scatterplot=0, bins=100,
                    variable=:exploit_tot, alpha=1, markersize=2, comparisonline=true,
                    scale=1)
    df = df[df[!,:capac].>=mincapac, :]
    plotly()
    blank = "<span style='color:white; font-size: 1px;'>q</span>"
    if scatterplot == 0
        p = histogram(df[!,variable], bins=range(0, 25.5, length=bins),
            xlabel="Exploited area per municipality/county [%]",
            ylabel="Number of municipalities/counties<br>$blank",
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            size=scale.*(800,550), legend=false)  # title=variable
        comparisonline && plot!(100*area_mult*[1, 1], collect(ylims(p)), line=(3, :dash))
        display(p)
        return df
    end
    # Filter out NaNs so scatter plots don't get messed up.
    # https://github.com/JuliaPlots/Plots.jl/issues/3258
    df = df[.!(isnan.(df[!,:capac]) .| isnan.(df[!,:class]) .| isnan.(df[!,variable])), :]
    if scatterplot == 1
        p = scatter(df[!,variable], df[!,:capac], 
            xlabel="Exploited area per municipality/county [%]",
            ylabel="MW",
            size=scale.*(800,550), legend=false, colorbar=true,
            markersize=scale*df[!,:class].^(1+markersize/10),
            marker_z=df[!,:masked], color=:watermelon, alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable]*area_mult,digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        display(p)
    elseif scatterplot == 2
        p = scatter(df[!,variable], df[!,:class], 
            xlabel="Exploited area per municipality/county [%]", xlims=(-0.5,25.5),
            ylabel="Mean wind class",
            size=scale.*(800,550), legend=false, colorbar=true,
            markersize=scale*df[!,:capac].^(markersize/7),
            markerstrokecolor=RGBA(0,0,0,.2), markerstrokewidth=1,
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            marker_z=df[!,:masked], color=:plasma, alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable]*area_mult,digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        display(p)
    elseif scatterplot == 3
        p = scatter(df[!,variable], df[!,:class], 
            xlabel="Exploited area per municipality/county [%]", xlims=(-0.5,25.5),
            ylabel="Mean wind class",
            size=scale.*(800,550), legend=false, colorbar=true,
            markersize=markersize*scale,
            markerstrokecolor=RGBA(0,0,0,.2), markerstrokewidth=0,
            tickfont=14*scale, guidefont=14*scale, left_margin=30px,
            color=RGB(0,0,0), alpha=alpha,
            hover=df[!,:region].*"<br>".*string.(df[!,:capac]).*" MW".*
            "<br>exploited = ".*string.(round.(df[!,variable]*area_mult,digits=1)).*"%".*
            "<br>class = ".*string.(df[!,:class]).*"<br>masked = ".*string.(df[!,:masked]))
        display(p)
    end 
    return df
end

function minimumpositive(arr::AbstractArray)
    min = typemax(eltype(arr))
    index = 0
    for (i, elem) in enumerate(arr)
        if elem > 0 && elem < min
            min = elem
            index = i
        end
    end
    return index, min
end

function closestturbine(region_abbrev)
    df = DataFrame!(CSV.File(in_datafolder("turbines_$region_abbrev.csv")))
    len = size(df, 1)
    m_per_degree = pi*6371*2/360*1000
    closest_dist = zeros(len)
    rotordiams = zeros(Int, len)
    for i = 1:len
        row = df[i, :]
        _, min = minimumpositive((df[:lat] .- row[:lat]).^2 + (df[:lon] .- row[:lon]).^2)
        closest_dist[i] = m_per_degree * sqrt(min)
        rotordiams[i] = row[:rotordiam]
    end
    # display(histogram(closest_dist[closest_dist .< 3000], bins=100))
    rotordist = closest_dist./rotordiams
    display(histogram(rotordist[rotordist .< 25], bins=1000))
    df, closest_dist, rotordiams
end

end # module
