module WindSites

using Proj4, XLSX, CSV, DataFrames, Dates, Pkg.TOML, Plots, MAT,
        Loess, SmoothingSplines, Plots.PlotMeasures, CategoricalArrays, PrettyTables

export openmap, readusa, readdk, readuk, readse, readde, readturbinedata, shapefile2csv,
    scatterplots_model, plotdist, fill_missing_rotordiams!, closestturbine, turbines2parks,
    analyze_protected, analyze_natura2000, analyze_landtype, grouppopulationdensity, groupwindspeeds,
    readfarms, savefarms, read_irena, stats

include("regional_level_GIS.jl")
include("turbine_level_GIS.jl")

function openmap(df::DataFrame, turbinenumber::Int)
    openmap(df, turbinenumber, :google)
    openmap(df, turbinenumber, :bing)
end
    
function openmap(df::DataFrame, turbinenumber::Int, source::Symbol)
    lon, lat = df[turbinenumber, :lon], df[turbinenumber, :lat]
    if source == :google
        # url = "https://www.google.com/maps/@?api=1\"&\"map_action=map\"&\"basemap=satellite\"&\"center=$lat%2C$lon\"&\"zoom=18"   # satellite map but no pin
        url = "https://www.google.com/maps/search/?api=1\"&\"query=$lat%2C$lon\"&\"basemap=satellite\"&\"zoom=18"   # map with pin but not satellite
        # url = "http://maps.google.com/maps?t=k\"&\"q=loc:$lat+$lon"   # used to do both but no longer places pin
        # https://stackoverflow.com/questions/47038116/google-maps-url-with-pushpin-and-satellite-basemap
        # https://stackoverflow.com/questions/60219254/show-location-marker-in-new-browser-window
        c = Cmd(`cmd /c start \"\" $url`, windows_verbatim=true)
        run(c)
        url = "https://www.google.com/maps/@?api=1\"&\"map_action=map\"&\"basemap=satellite\"&\"center=$lat%2C$lon\"&\"zoom=18" 
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

function readall()
    cols = [:lon, :lat, :capac, :year, :onshore, :elec2018]
    df_DK = DataFrame(CSV.File(in_datafolder("turbines_DK.csv")))[:, cols]
    df_SE = DataFrame(CSV.File(in_datafolder("turbines_SE.csv")))[:, cols[1:5]]
    df_DE = DataFrame(CSV.File(in_datafolder("turbines_DE.csv")))[:, cols[1:4]]
    df_USA = DataFrame(CSV.File(in_datafolder("turbines_USA.csv")))[:, cols[1:4]]
    df_SE[:, :elec2018] .= missing    # MWh/year
    df_DE[:, :elec2018] .= missing    # MWh/year
    df_DE[:, :onshore] .= missing    # MWh/year
    df_USA[:, :elec2018] .= missing    # MWh/year
    df_USA[:, :onshore] .= missing    # MWh/year
    df_SE[!,:year] = [ismissing(d) ? 2000 : year(d) for d in df_SE[!,:year]]
    df_DK.country .= "Denmark"
    df_SE.country .= "Sweden"
    df_DE.country .= "Germany"
    df_USA.country .= "USA"
    return vcat(df_DK, df_SE, df_DE, df_USA)
end

function readfarms(winddir = "C:/Users/niclas/Downloads/wind data")
    df = DataFrame(CSV.File("$winddir/Windfarms_World_20220407_fixed.csv"; quotechar='\'', missingstring=["#ND", ""]))
    # ["ID (#ND = no data)", "Continent", "ISO code (Code ISO 3166.1)", "Country", "State code", "Area", "City", "Name", "2nd name",
    #     "Latitude (WGS84)", "Longitude (WGS84)", "Altitude/Depth (m)", "Location accuracy (Yes = accurate location)", "Offshore - Shore distance (km)",
    #     "Manufacturer", "Turbine", "Hub height (m)", "Number of turbines", "Total power (kW)", "Developer", "Operator", "Owner",
    #     "Commissioning date (Format: yyyy or yyyymm)", "Status", "Decommissioning date (Format: yyyy or yyyymm)", "Link", "Update"]
    newnames = [:id, :continent, :iso, :country, :state, :area, :city, :name, :name2, :lat, :lon, :altitude, :accurate_location, :offshore,
        :manufacturer, :turbine, :hubheight, :num_turbines, :capac, :developer, :operator, :owner, :startdate, :status, :enddate, :link, :updated]
    rename!(df, newnames)
    df[!, [:lat, :lon]] = [ismissing(x) ? missing : parse(Float64, replace(x, "," => ".")) for x in Array(df[!, [:lat, :lon]])]
    df.year = [ismissing(x) ? missing : parse(Int, x[1:4]) for x in df.startdate]
    df.endyear = [ismissing(x) ? missing : parse(Int, x[1:4]) for x in df.enddate]
    df.month = [ismissing(x) ? missing : (m = match(r".*/(\d+)", x)) === nothing ? missing : parse(Int, m[1]) for x in df.startdate]
    df.endmonth = [ismissing(x) ? missing : (m = match(r".*/(\d+)", x)) === nothing ? missing : parse(Int, m[1]) for x in df.enddate]
    df.accurate_location = [x == "Yes" ? true : x == "No" ? false : missing for x in df.accurate_location]
    df.onshore = [x == "No" ? true : x == "Yes" ? false : missing for x in df.offshore]
    df.altitude = [ismissing(x) ? missing : (m = match(r"(\d+)/(\d+)", x)) === nothing ? round(Int, parse(Float64, x)) :
                            round(Int, mean(parse.(Int, [m[1], m[2]]))) for x in df.altitude]
    replacecountries = ["United-Kingdom" => "United Kingdom", "New-Zealand" => "New Zealand"]
    replace!(df.country, replacecountries...)
    filter!(row -> row.status == "Production", df)
    select!(df, [newnames[2:13]; :onshore; newnames[15:22]; :year; :month; :status; :endyear; :endmonth])
    return df
end

function savefarms(winddir = "C:/Users/niclas/Downloads/wind data")
    wf = readfarms(winddir)
    CSV.write("$winddir/windfarms.csv", wf)
end

function read_irena(winddir = "C:/Users/niclas/Downloads/wind data")
    df = DataFrame(CSV.File("$winddir/IRENA ELECCAP_20220411-085642.csv", header=3, missingstring=".."))
    # "Country/area", "Technology", "Grid connection", "Year", "Installed electricity capacity by country/area (MW)"
    newnames = [:country, :onshore, :gridconnected, :year, :capac]
    rename!(df, newnames)
    df.onshore = (df.onshore .== "Onshore wind energy")
    df.gridconnected = (df.gridconnected .== "On-grid")
    replacecountries = ["United Kingdom of Great Britain and Northern Ireland" => "United Kingdom", "United States of America" => "USA"]
    replace!(df.country, replacecountries...)
    # rows = df.gridconnected .&& df.year .== 2020
    # onshore2020 = Dict(row.country => row.capac for row in eachrow(df[rows .&& df.onshore, :]) if !ismissing(row.capac))
    # offshore2020 = Dict(row.country => row.capac for row in eachrow(df[rows .&& .!df.onshore, :]) if !ismissing(row.capac))
    return df
end

function stats(onoff)
    countries = ["Denmark", "Sweden", "Ireland", "Germany", "Norway", "Spain", "Portugal", "Finland", "Uruguay", "Belgium", "Greece",
                    "Netherlands", "United Kingdom", "Australia", "Austria", "USA", "Canada", "Luxembourg", "France", "Estonia"]
    wf = readfarms()
    dff = readall()
    dfi = read_irena()
    headings = ["country", "WF_2020", "IRENA_2020", "WF_share", "FRT_2020", "FRT_share", "WF_has_year", "WF_has_location", "WF_accurate_location"]
    data = reshape([], 0, length(headings))
    for country in countries
        onshore = onoff == "onshore" ? dfi.onshore : onoff == "offshore" ? .!dfi.onshore : ones(Bool, size(dfi, 1)) 
        capac_irena = sum(skipmissing(dfi[dfi.country .== country .&& onshore .&& dfi.gridconnected .&& dfi.year .== 2020, :capac])) / 1000 |> x -> round(x, digits=2)
        if country in ["Denmark", "Sweden", "Germany", "USA"]
            onshore = onoff == "onshore" ? dff.onshore : onoff == "offshore" ? .!dff.onshore : ones(Bool, size(dff, 1)) 
            capac_frt = sum(dff[dff.country .== country .&& .!ismissing.(onshore) .&& onshore .&& .!ismissing.(dff.year) .&& dff.year .<= 2020, :capac]) / 1e6 |> x -> round(x, digits=2)
            frt_share = round(100*capac_frt/capac_irena, digits=1)
        else
            capac_frt, frt_share = ".", "."
        end
        onshore = onoff == "onshore" ? wf.onshore : onoff == "offshore" ? .!wf.onshore : ones(Bool, size(wf, 1)) 
        df = wf[wf.country .== country .&& .!ismissing.(onshore) .&& onshore .&& .!ismissing.(wf.capac), :]
        capac_wf = sum(df[.!ismissing.(df.year) .&& df.year .<= 2020, :capac]) / 1e6 |> x -> round(x, digits=2)
        wf_share = round(100*capac_wf/capac_irena, digits=1)
        capac_has_year = sum(df[.!ismissing.(df.year), :capac])
        capac_has_location = sum(df[.!ismissing.(df.lat), :capac])
        capac_accurate_location = sum(df[df.accurate_location, :capac])
        capac_tot = sum(df[:, :capac])
        has_year = round(100*capac_has_year/capac_tot, digits=1)
        has_location = round(100*capac_has_location/capac_tot, digits=1)
        accurate_location = round(100*capac_accurate_location/capac_tot, digits=1)
        if capac_irena == 0
            wf_share = "."
        end
        if capac_tot == 0
            has_year, has_location, accurate_location = ".", ".", "."
        end
        data = [data; country capac_wf capac_irena wf_share capac_frt frt_share has_year has_location accurate_location]
    end
    pretty_table(data, header=(headings, ["", "GW", "GW", "% (irena)", "GW", "% (irena)", "% (capac)", "% (capac)", "% (capac)"]))
end

function readusa()
    # US Wind Turbine Database
    # https://eerscmap.usgs.gov/uswtdb/
    datadir = joinpath(@__DIR__, "..", "data")
    df = DataFrame(CSV.File("$datadir/USA_uswtdb_v3_1_20200717.csv"))
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
    df = DataFrame(CSV.File("$datadir/UK_REPD-June-2020-update.csv"))[:, cols]
    rename!(df, [:type, :capac_park, :capac, :nturbines, :height, :status, :lon, :lat, :year])
    delete!(df, .!startswith.(df[!, :type], "Wind"))
    delete!(df, df[:, :status] .!= "Operational")
    select!(df, Not(:status))
    df[!, :type] = replace.(df[!, :type], "Wind " => "")
    df[!, :nturbines] = parse.(Int, df[!, :nturbines])
    df[!, :height] = parse_not_missing.(Int, df[!, :height])
    df[!, :capac] = parse.(Float64, df[!, :capac])
    df[!, :capac_park] = parse.(Float64, replace.(df[!, :capac_park], ", " => ""))
    df[!, :lon] = parse.(Int, replace.(df[!, :lon], ", " => ""))
    df[!, :lat] = parse.(Int, replace.(df[!, :lat], ", " => ""))
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
    df = DataFrame(XLSX.readtable("$datadir/DK_anlaegprodtilnettet_0.xlsx", "IkkeAfmeldte-Existing turbines",
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
    df = DataFrame(XLSX.readtable("$datadir/SE_Vindbrukskollen_export_allman_Prod.xlsx", "Vindkraftverk",
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
    df = DataFrame(CSV.File("$datadir/DE_windpower.csv"))[!, 2:7]
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
    df = DataFrame(CSV.File(in_datafolder("turbines_$region_abbrev.csv")))
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
