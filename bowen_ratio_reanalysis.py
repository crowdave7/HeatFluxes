"""Import necessary modules for this code."""
import cartopy as cart
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris
import iris.analysis
import iris.coord_categorisation
import iris.plot as iplt
import matplotlib
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np
import os


def reanalysis_bowen(reanalysis_type, season_name):
    """Take the input variables, and find the paths to the relevant regridded sensible and latent heat flux model files."""
    """Plot these up as a spatial map for the given season."""
    """Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles"

    """Find the paths to the files containing the model data"""
    model_file_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            for j in reanalysis_type:
                if j in path and ('hfss' in path or 'hfls' in path):
                    model_file_paths = np.append(model_file_paths, path)

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    print model_file_paths

    """Load the sensible and latent heat flux data from each model into a cubelist."""

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    """Load the two cubes."""
    cubes = iris.load(model_file_paths)
    count = 0

    """Constrain the time range of the cubes."""
    for i in cubes:
        with iris.FUTURE.context(cell_datetime_objects=True):
            cubes[count] = i.extract(time_range)
            time_points = cubes[count].coord('time').points
            times = cubes[count].coord('time').units.num2date(time_points)
            model_id = reanalysis_type
            print model_id
            print len(times)
            print times[0]
            print times[-1]
            count +=1

    """Take the mean over time for each cube in the cubelist."""
    cube_id = 0

    """For each cube (for the sensible and latent heat flux data for one model)"""
    for i in cubes:

        """Slice the regridded cube down to the African domain."""
        cube = i.intersection(latitude=(-40, 40), longitude=(-30, 70))

        time_points = cube.coord('time').points
        times = cube.coord('time').units.num2date(time_points)

        """Plot up a map for the file."""

        """ If the input month is defined as the whole year,"""
        if season_name == 'Climatology':

            """Take the mean over the cube."""
            reanalysis_data = cube.collapsed('time', iris.analysis.MEAN)

            """Replace the cube in the cubelist with the new cube."""
            cubes[cube_id] = reanalysis_data

            """Add 1 to the count variable."""
            cube_id +=1

        """ If the input month is defined as a season,"""
        if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

            """Create two coordinates on the cube to represent seasons and the season year."""
            iris.coord_categorisation.add_season(i, 'time', name='clim_season')
            iris.coord_categorisation.add_season_year(i, 'time', name='season_year')

            """Aggregate the data by season and season year."""
            seasonal_means = i.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

            """Constrain the data by the input season. Decapitalise the input variable."""
            constraint = iris.Constraint(clim_season=season_name.lower())
            reanalysis_data_season = seasonal_means.extract(constraint)

            """Take the mean over the cube."""
            reanalysis_data = reanalysis_data_season.collapsed('time', iris.analysis.MEAN)

            """Replace the cube in the cubelist with the new cube."""
            cubes[cube_id] = reanalysis_data

            """Add 1 to the count variable."""
            cube_id +=1

    """Select the numerator (latent heat flux)."""
    variable_name_0 = str(cubes[0].long_name)
    variable_name_1 = str(cubes[1].long_name)
    if "Latent" in variable_name_0 or "latent" in variable_name_0:
        numerator = cubes[0]
    elif "Latent" in variable_name_1 or "latent" in variable_name_1:
        numerator = cubes[1]

    """Calculate the denominator (latent + sensible heat flux)."""
    denominator = iris.analysis.maths.add(cubes[0], cubes[1])

    """Calculate the evaporative fraction (latent/(latent+sensible))."""
    evap_fraction = iris.analysis.maths.divide(numerator, denominator)

    """Plot the figure."""
    fig = plt.figure()

    """Import coastlines and lake borders. Set the scale to 10m, 50m or 110m resolution for more detail."""
    coastline = cart.feature.NaturalEarthFeature(category='physical', name='coastline', scale='110m', facecolor='none')
    lake_borders = cart.feature.NaturalEarthFeature(category='physical', name='lakes', scale='110m', facecolor='none')

    """Import country borders."""
    shapefile = shapereader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries')
    reader = shapereader.Reader(shapefile)
    country_borders = reader.records()

    """Remove iris warning message."""
    iris.FUTURE.netcdf_promote = True

    """Define the contour levels for the input variables."""
    contour_levels = np.arange(0.4, 1.05, 0.05)

    """Define the colour map and the projection."""
    cmap = matplotlib.cm.get_cmap('YlGnBu')
    crs_latlon = ccrs.PlateCarree()

    """Plot the map using cartopy, and add map features."""
    ax = plt.subplot(111, projection=crs_latlon)
    ax.set_extent([-22, 62, -22, 12], crs=crs_latlon)
    contour_plot = iplt.contourf(evap_fraction, contour_levels, cmap=cmap, extend='both')
    ax.add_feature(coastline, zorder=5, edgecolor='k', linewidth=2)
    ax.add_feature(lake_borders, zorder=5, edgecolor='k', linewidth=1)
    for i in country_borders:
        ax.add_geometries(i.geometry, ccrs.PlateCarree(), edgecolor="black", facecolor="None")

    """Define gridlines."""
    gridlines = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='black', linewidth=0.4, linestyle='--')
    gridlines.xlabels_top = False
    gridlines.xlabels_bottom = True
    gridlines.ylabels_left = True
    gridlines.ylabels_right = False
    gridlines.xlines = True
    gridlines.ylines = True
    gridlines.xformatter = LONGITUDE_FORMATTER
    gridlines.yformatter = LATITUDE_FORMATTER
    gridlines.xlocator = mticker.FixedLocator(np.arange(-40, 100, 20))
    gridlines.ylocator = mticker.FixedLocator(np.arange(-50, 70, 10))

    """Add a colour bar, with ticks and labels."""
    colour_bar = plt.colorbar(contour_plot, orientation='horizontal', pad=0.1, aspect=40)
    colour_bar.set_ticks(np.arange(0.4, 1.05, 0.05))
    colour_bar.set_ticklabels(np.arange(0.4, 1.05, 0.05))
    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)

    variable_name = "Evaporative Fraction"
    colour_bar.set_label((variable_name), fontsize=10)

    if reanalysis_type == ["erai"]:
        title = "ERA-Interim 1979-2008 "

    if reanalysis_type == ["merra2"]:
        title = "MERRA-2 1979-2008"

    if reanalysis_type == ["doe"]:
        title = "NCEP DOE-2 1979-2008 "

    """Add a title."""
    plt.title(title+""+season_name+"", fontsize=10)

    """Save the figure, close the plot and print an end statement."""
    fig.savefig("Evap_Frac_"+season_name+"_"+reanalysis_type[0]+".png")
    plt.close()
    print "Reanalysis Plot done"


#reanalysis_bowen(["erai"], "SON")
