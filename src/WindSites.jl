module WindSites

using Proj4, XLSX, CSV, DataFrames, Dates, GDAL_jll, Pkg.TOML, Plots

export openmap, readusa, readdk, readuk, readse, readturbinedata, shapefile2csv, scatterplots

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

function readse()
    # Länsstyrelsen: Vindbrukskollen
    # https://vbk.lansstyrelsen.se/  (click "Excel-export")
    # shapefile: https://ext-dokument.lansstyrelsen.se/gemensamt/geodata/ShapeExport/LST.VKOLLEN_VINDKRAFTVERK.zip 
    # shapefile2csv("LST.VKOLLEN_VINDKRAFTVERK/LST_VKOLLEN_VVERK.shp")

    datadir = joinpath(@__DIR__, "..", "data")
    df = DataFrame(CSV.File("$datadir/SE_Vindbrukskollen_export_allman_Prod.csv"))
    select!(df, ["Status", "Placering", "E-Koordinat", "N-Koordinat", "Totalhöjd (m)", "Navhöjd (m)",
        "Rotordiameter (m)", "Maxeffekt (MW)", "Uppfört", "Fabrikat", "Modell"])
    rename!(df, [:status, :type, :lon, :lat, :height, :hubheight, :rotordiam, :capac, :year, :brand, :model])
    delete!(df, df[:, :status] .!= "Uppfört")
    select!(df, Not(:status))
    df[!, :onshore] = startswith.(df[!, :type], "Land")     # all turbines are marked "Land" or "Vatten"
    select!(df, Not(:type))
    delete!(df, ismissing.(df[:, :capac]))
    delete!(df, df[:, :capac] .== 0)
    df[!, :capac] = Int.(df[!, :capac] * 1000)
    df[!, :year] = [d == "Jan 01, 1900" ? missing : Date(d, dateformat"u d, y") for d in df[!, :year]]
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

function readturbinedata()
    df_dk = readdk()
    df_usa = readusa()
    df_uk = readuk()
    df_se = readse()
    df_de = readde()

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

function scatterplots(gisregion, type=[:total])
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

    df = innerjoin(dfg, dfm, on=:region)
    select!(df, [:region, :capac, :offcapac, :elec2018, :mcap, :moffcap, :melec, :moffelec, :mclass, :moffclass])
    plotly()
    println(type)
    if :onshore in type
        s = scatter(df[:capac], df[:mcap], markersize=df[:mclass]*2, title="Onshore capacity",
            hover=df[:region].*" class ".*string.(df[:mclass]),
            xlabel="GW (real)", ylabel="GW (model)", size=(800,650))
        display(s)
    elseif :offshore in type
        s = scatter(df[:offcapac], df[:moffcap], markersize=df[:moffclass]*2, title="Offshore capacity",
            hover=df[:region].*" class ".*string.(df[:moffclass]),
            xlabel="GW (real)", ylabel="GW (model)", size=(800,650))
        display(s)
    elseif :total in type
        s = scatter(df[:capac]+df[:offcapac], df[:mcap]+df[:moffcap], markersize=df[:mclass]*2, title="Total capacity",
            hover=df[:region].*" class ".*string.(df[:mclass]),
            xlabel="GW (real)", ylabel="GW (model)", size=(800,650))
        display(s)
    elseif :elec in type
        s = scatter(df[:elec2018], df[:melec], markersize=df[:mclass], xlabel="GWh (real)", ylabel="GWh (model)")
        display(s)
    end
    df
 end


end # module
