function analyze_protected(; firstyear=1978, lastyear=2021)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame

    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear)
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    df.countryname = countries[df.country]
    df.onshore = df.landtype .> 0
    protected_names = [
       "(Ia) Strict Nature Reserve"
       "(Ib) Wilderness Area"
       "(II) National Park"
       "(III) Natural Monument"
       "(IV) Habitat/Species Management"
       "(V) Protected Landscape/Seascape"
       "(VI) Managed Resource Protected Area"
       "Not Reported"
       "Not Applicable"
       "Not Assigned"
    ]
    protected_type = Dict(i => n for (i,n) in enumerate(protected_names))
    protected_type[0] = "UNPROTECTED"

    gdf = groupby(df, [:countryname, :onshore, :protected])
    cdf = combine(gdf,
            :area => (a -> sum(a)/1000) => :area,
            [:turbinecapac, :inyears] => ((c,y) -> sum(c.*y)/1000) => :capac
        )
    gdf_tot = groupby(cdf, [:countryname, :onshore])
    cdf_tot = combine(gdf_tot, :area => sum, :capac => sum)
    
    cdf.turbinedensity = cdf.capac ./ cdf.area
    cdf.protected_type = getindex.(Ref(protected_type), cdf.protected)
    cdf_tot.turbinedensity = cdf_tot.capac_sum ./ cdf_tot.area_sum

    sort!(cdf, [order(:countryname, by = x -> findfirst(countries.==x)), :protected, order(:onshore, rev=true)])
    sort!(cdf_tot, [order(:countryname, by = x -> findfirst(countries.==x)), order(:onshore, rev=true)])
    # CSV.write("protectedinfo.csv", cdf)
    # CSV.write("protectedinfo_tot.csv", cdf_tot)

    filename = "protected WDPA.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [cdf[:,i] for i=1:size(cdf,2)], names(cdf), anchor_cell=XLSX.CellRef("A4"))
        XLSX.writetable!(sheet, [cdf_tot[:,i] for i=1:size(cdf_tot,2)], names(cdf_tot), anchor_cell=XLSX.CellRef("L4"))
    end

    return cdf, cdf_tot
end

function analyze_natura2000(; firstyear=1978, lastyear=2021)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame

    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear)
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    df.countryname = countries[df.country]
    df.onshore = df.landtype .> 0
    protected_names = [
       "A: SPAs"
       "B: SCIs and SACs"
       "C: both categories A + B"
    ]
    protected_type = Dict(i => n for (i,n) in enumerate(protected_names))
    protected_type[0] = "UNPROTECTED"

    gdf = groupby(df, [:countryname, :onshore, :natura2000])
    cdf = combine(gdf,
            :area => (a -> sum(a)/1000) => :area,
            [:turbinecapac, :inyears] => ((c,y) -> sum(c.*y)/1000) => :capac
        )
    gdf_tot = groupby(cdf, [:countryname, :onshore])
    cdf_tot = combine(gdf_tot, :area => sum, :capac => sum)
    
    cdf.turbinedensity = cdf.capac ./ cdf.area
    cdf.protected_type = getindex.(Ref(protected_type), cdf.natura2000)
    cdf_tot.turbinedensity = cdf_tot.capac_sum ./ cdf_tot.area_sum

    sort!(cdf, [order(:countryname, by = x -> findfirst(countries.==x)), :natura2000, order(:onshore, rev=true)])
    sort!(cdf_tot, [order(:countryname, by = x -> findfirst(countries.==x)), order(:onshore, rev=true)])

    filename = "protected Natura2000.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [cdf[:,i] for i=1:size(cdf,2)], names(cdf), anchor_cell=XLSX.CellRef("A4"))
        XLSX.writetable!(sheet, [cdf_tot[:,i] for i=1:size(cdf_tot,2)], names(cdf_tot), anchor_cell=XLSX.CellRef("L4"))
    end
    # CSV.write("protectedinfo.csv", cdf)
    # CSV.write("protectedinfo_tot.csv", cdf_tot)

    return cdf, cdf_tot
end

function analyze_landtype(; firstyear=1978, lastyear=2021, minspeed=0)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame

    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear) .& (df.windspeed .>= minspeed)
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    df.countryname = countries[df.country]
    landtype_names = [
        "Water" 
        "Evergreen Needleleaf Forests"
        "Evergreen Broadleaf Forests"
        "Deciduous Needleleaf Forests"
        "Deciduous Broadleaf Forests"
        "Mixed Forests"
        "Closed Shrublands"
        "Open Shrublands"
        "Woody Savannas"
        "Savannas"
        "Grasslands"
        "Permanent Wetlands"
        "Croplands"
        "Urban"
        "Cropland/Natural"
        "Snow/Ice"
        "Barren"
    ]
    landtype = Dict(i-1 => n for (i,n) in enumerate(landtype_names))

    gdf = groupby(df, [:countryname, :landtype])
    cdf = combine(gdf,
            :area => (a -> sum(a)/1000) => :area,
            [:turbinecapac, :inyears] => ((c,y) -> sum(c.*y)/1000) => :capac
        )
    gdf_tot = groupby(cdf, :countryname)
    cdf_tot = combine(gdf_tot, :area => sum, :capac => sum)
    cdf.turbinedensity = cdf.capac ./ cdf.area
    cdf.landtype_name = getindex.(Ref(landtype), cdf.landtype)
    cdf_tot.turbinedensity = cdf_tot.capac_sum ./ cdf_tot.area_sum

    sort!(cdf, [order(:countryname, by = x -> findfirst(countries.==x)), :landtype])
    sort!(cdf_tot, [order(:countryname, by = x -> findfirst(countries.==x))])

    filename = "landtypes minspeed $minspeed.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [cdf[:,i] for i=1:size(cdf,2)], names(cdf), anchor_cell=XLSX.CellRef("A4"))
        XLSX.writetable!(sheet, [cdf_tot[:,i] for i=1:size(cdf_tot,2)], names(cdf_tot), anchor_cell=XLSX.CellRef("L4"))
    end

    return cdf, cdf_tot
end

function grouppopulationdensity(country; firstyear=1978, lastyear=2021)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame
    df = df[df.country .== country, :]
    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear)
    df.nturbines = df.nturbines .* df.inyears
    df.capac = df.turbinecapac/1000 .* df.inyears
    pops = [[i*10.0^j for j = -1:3 for i in [1,2,5]]; 100_000]
    popmin = [0; pops[1:end-1]]
    df.poprange = CategoricalArrays.cut(df.popdens, pops, extend=true)
    gdf = groupby(df, :poprange)
    cdf = combine(gdf, :nturbines .=> sum, :capac .=> sum, nrow => :pixels, :area .=> sum, 
        renamecols = false)    
    allranges = CategoricalArrays.cut(popmin, pops, extend=true)
    df_out = similar(cdf, length(allranges))
    df_out.poprange = allranges
    df_out[:,2:end] .= 0
    insertcols!(df_out, 2, :popmin => popmin, :popmax => pops)
    indexes = [findfirst(allranges .== r) for r in cdf.poprange]
    df_out[indexes, 4:end] = cdf[:, 2:end]
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    
    filename = "wind_popdens 1978-2021.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$(countries[country]) $firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [df_out[:,i] for i=2:size(df_out,2)], names(df_out)[2:end], anchor_cell=XLSX.CellRef("A1"))
    end
    # CSV.write(filename, df_out[:, 2:end])
end

# just for the histogram data Fredrik wanted
function groupwindspeeds(country; firstyear=1978, lastyear=2021, usemiuu=false)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame
    df = df[df.country .== country, :]
    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear)
    df.is_onshore = (df.landtype .> 0)
    df.is_offshore = (df.landtype .== 0)
    df.nturbines = df.nturbines .* df.inyears
    df.nturbines_on = df.nturbines .* df.is_onshore
    df.nturbines_off = df.nturbines .* df.is_offshore
    df.capac = df.turbinecapac/1000 .* df.inyears
    df.capac_on = df.capac .* df.is_onshore
    df.capac_off = df.capac .* df.is_offshore
    df.area_on = df.area .* df.is_onshore
    df.area_off = df.area .* df.is_offshore
    df.speed = usemiuu ? df.miuu : df.windspeed
    df.windspeed_range = CategoricalArrays.cut(df.speed, 0:.25:21, extend=true)
    gdf = groupby(df, :windspeed_range)
    cdf = combine(gdf, [:nturbines, :nturbines_on, :nturbines_off] .=> sum, 
        [:capac, :capac_on, :capac_off] .=> sum,
        nrow => :pixels, [:is_onshore, :is_offshore] .=> sum .=> [:pixels_on, :pixels_off],
        [:area, :area_on, :area_off] .=> sum, 
        renamecols = false)
    allranges = CategoricalArrays.cut(0:.25:20.75, 0:.25:21, extend=true)
    df_out = similar(cdf, length(allranges))
    df_out.windspeed_range = allranges
    df_out[:,2:end] .= 0
    indexes = [findfirst(allranges .== r) for r in cdf.windspeed_range]
    df_out[indexes, :] = cdf
    insertcols!(df_out, 2, :windspeed_low => 0:.25:20.75, :windspeed_high => .25:.25:21)

    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]

    filename = "winddata 1978-2021.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$(countries[country]) $(usemiuu ? "MIUU " : "")$firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [df_out[:,i] for i=2:size(df_out,2)], names(df_out)[2:end], anchor_cell=XLSX.CellRef("A1"))
    end
    # CSV.write(filename, df_out)
end
