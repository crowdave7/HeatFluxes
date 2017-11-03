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


def plot_bowen(list_of_models, model_type, season_name):
    """Take the input variables, and find the paths to the relevant regridded sensible and latent heat flux model files."""
    """Plot these up as a spatial map for the given season."""
    """Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles"

    """Find the paths to the files containing the model data"""
    model_file_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            #print path
            for j in list_of_models:
                for char in '/':
                    j = j.replace(char,'')
                if j in path and model_type in path and ('hfss' in path or 'hfls' in path):
                    model_file_paths = np.append(model_file_paths, path)

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    print model_file_paths

    """Load the sensible and latent heat flux data from each model into a cubelist."""

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    """Load the data from the model file paths into a cube. Constrain the input years"""
    """Print the model ID, length of time dimension, and first and last model dates."""

    """Create a blank list to contain all the model data for the evaporative fraction."""
    cubes_ef = []

    """ For each model compute the evaporative fraction in an iris cube and plot."""
    for j in list_of_models:
        for char in '/':
            j = j.replace(char,'')

        """Find the relevant model paths, and load into a cubelist."""
        paths_for_this_model = [k for k in model_file_paths if j in k]
        cubes = iris.load(paths_for_this_model, ['surface_upward_latent_heat_flux', 'surface_upward_sensible_heat_flux'])

        """Constrain the years. Print the model ID, length of time dimension, and first and last model dates."""
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

        """Take the mean over time for each cube in the cubelist."""
        cube_id = 0

        """Print the model id from one of the cubes."""
        model_id = cubes[0].attributes['model_id']

        """For each cube,"""
        for regridded_model_data in cubes:

            """ If the input month is defined as the whole year,"""
            if season_name == 'Climatology':

                """Take the mean over the cube."""
                model_data = regridded_model_data.collapsed('time', iris.analysis.MEAN)

                """Replace the cube in the cubelist with the new cube."""
                cubes[cube_id] = model_data

                """Add 1 to the count variable."""
                cube_id +=1

            """ If the input month is defined as a season,"""
            if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

                """Create two coordinates on the cube to represent seasons and the season year."""
                iris.coord_categorisation.add_season(regridded_model_data, 'time', name='clim_season')
                iris.coord_categorisation.add_season_year(regridded_model_data, 'time', name='season_year')

                """Aggregate the data by season and season year."""
                seasonal_means = regridded_model_data.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

                """Constrain the data by the input season. Decapitalise the input variable."""
                constraint = iris.Constraint(clim_season=season_name.lower())
                regridded_model_data_season = seasonal_means.extract(constraint)

                """Take the mean over the cube."""
                model_data = regridded_model_data_season.collapsed('time', iris.analysis.MEAN)

                """Replace the cube in the cubelist with the new cube."""
                cubes[cube_id] = model_data

                """Add 1 to the count variable."""
                cube_id +=1

        """Now compute the evaporative fraction."""

        """Select the numerator (latent heat flux)."""
        variable_name_0 = str(cubes[0].standard_name)
        variable_name_1 = str(cubes[1].standard_name)
        if "latent" in variable_name_0:
            numerator = cubes[0]
        if "latent" in variable_name_1:
            numerator = cubes[1]

        """Calculate the denominator (latent + sensible heat flux)."""
        denominator = iris.analysis.maths.add(cubes[0], cubes[1])

        """Calculate the evaporative fraction (latent/(latent+sensible))."""
        evap_fraction = iris.analysis.maths.divide(numerator, denominator)

        cubes_ef = np.append(cubes_ef, evap_fraction)

    ensemble_mean = sum(cubes_ef) / float(len(cubes_ef))

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
    cmap = matplotlib.cm.get_cmap('coolwarm_r')
    crs_latlon = ccrs.PlateCarree()

    """Plot the map using cartopy, and add map features."""
    ax = plt.subplot(111, projection=crs_latlon)
    ax.set_extent([-22, 62, -22, 12], crs=crs_latlon)
    contour_plot = iplt.contourf(ensemble_mean, contour_levels, cmap=cmap, extend='both')
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

    """Add a title."""
    plt.title('Ensemble Mean (AMIP) 1979-2008 ' ''+season_name+'', fontsize=10)

    """Save the figure, close the plot and print an end statement."""
    fig.savefig("Evap_Frac_"+season_name+"_EnsembleMean.png")
    plt.close()
    print "Ensemble Mean Plot done"


#plot_bowen(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "MIROC5", "MPI-ESM-MR", "MRI-CGCM3", "NorESM1-M/"], "amip", "SON")
