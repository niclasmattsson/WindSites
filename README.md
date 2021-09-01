# WindSites
 
Below is some example code to run functions in this package.

## Setup

* Create a folder to contain dataframe files exported from the GlobalEnergyGIS package and move those files inside. Call the folder GISdata or similar.
* Create a file `$HOMEDIR/.GlobalEnergyGIS_config` containing one line `datafolder = "D:/GISdata"` and replace that path with the path to the above folder.
* Optional: Set the system environment variable `JULIA_PKG_DEVDIR` if you want to make changes to this package.

## Read dataframes with turbine data

```
using WindSites   # and wait for everything to precompile 
df_dk, df_usa, df_uk, df_se, df_de = readturbinedata()
```

## Open Google and Bing maps to a turbine location

```
openmap(df_se, 20)   # show turbine number 20
```

## Municipal level histogram

Using `scatterplot=0` will plot a histogram.

```
df = plotdist(["USAGADM3"], area_mult=0.08, variable=:exploit_tot, mincapac=10, scatterplot=0, alpha=.15, markersize=5, line=:mean, bins=200, xmax=80)
```

## Municipal level scatter plot

Try `scatterplot=x` with x between 1 and 6 for different chart variants.

```
# scatterplot=1, four regions
df = plotdist(["USAGADM3", "GermanyGADM3", "SwedenGADM3", "Denmark83"], area_mult=0.08, variable=:exploit_tot, mincapac=10, scatterplot=1, alpha=.35, markersize=2, line=:mean, bins=100, xmax=30);

# scatterplot=2, four regions
df = plotdist(["USAGADM3", "GermanyGADM3", "SwedenGADM3", "Denmark83"], area_mult=0.08, variable=:exploit_tot, mincapac=10, scatterplot=2, alpha=.15, markersize=2, line=:mean, bins=100, xmax=30);

# scatterplot=3, four regions
df = plotdist(["USAGADM3", "GermanyGADM3", "SwedenGADM3", "Denmark83"], area_mult=0.08, variable=:exploit_tot, mincapac=10, scatterplot=3, alpha=.15, markersize=2, line=:mean, bins=100, xmax=30);

# scatterplot=4, four regions
df = plotdist(["USAGADM3", "GermanyGADM3", "SwedenGADM3", "Denmark83"], area_mult=0.08, variable=:exploit_tot, mincapac=10, scatterplot=4, alpha=.1, markersize=4, line=:mean, bins=100, xmax=30);

# scatterplot=5, one region
df = plotdist(["GermanyGADM3"], area_mult=0.08, variable=:exploit_tot, mincapac=10, scatterplot=5, alpha=.15, markersize=2, line=:mean, bins=100, xmax=30, legend=false)

# scatterplot=6, one region
df = plotdist(["USAGADM3", "GermanyGADM3", "SwedenGADM3", "Denmark83"], area_mult=0.08, variable=:exploit_tot, mincapac=10, scatterplot=6, alpha=.15, markersize=2, line=:mean, bins=100, xmax=30);
```

## Export WDPA protected area statistics

This will create an excel file in your Julia work folder. You can restrict turbines by year of installation.

```
cdf, cdf_tot = analyze_protected()

cdf, cdf_tot = analyze_protected(lastyear=2010)

cdf, cdf_tot = analyze_protected(firstyear=2011, lastyear=2020)
```

## Export Natura2000 protected area statistics

You can also restrict turbines by year of installation as above.

```
cdf, cdf_tot = analyze_natura2000()
```

## Export land type statistics

This will export a file of land types for turbines with an annual average wind speed above 5 m/s. You can also restrict turbines by year of installation as above.

```
cdf, cdf_tot = analyze_landtype(minspeed=5)
```

## Export population density statistics

This will export population density statistics for a given country/state index number (1-14). You can also restrict turbines by year of installation as above.

```
grouppopulationdensity(1)

grouppopulationdensity(2, firstyear=2011, lastyear=2020)

for r=1:14
    println(r)
    grouppopulationdensity(r)
    grouppopulationdensity(r, lastyear=2010)
    grouppopulationdensity(r, firstyear=2011)
end
```

## Export wind speed statistics

This will export wind speed statistics for a given country/state index number (1-14). You can also restrict turbines by year of installation as above. For Sweden (index=1), you can choose to use the MIUU atlas of wind speeds instead of the default Global Wind Atlas.

```
for r=1:14
    println(r)
    groupwindspeeds(r)
    groupwindspeeds(r, lastyear=2010)
    groupwindspeeds(r, firstyear=2011)
end

groupwindspeeds(1, usemiuu=true)
groupwindspeeds(1, lastyear=2010, usemiuu=true)
groupwindspeeds(1, firstyear=2011, usemiuu=true)
```
