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


def map(list_of_models, model_type, variable, season_name):
    """Take the input variables, and find the paths to the relevant model files."""
    """Plot these up as a spatial map for the given season."""
    """Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

    """If variable is pr, distinguish between pr and precipitable water to find model files."""
    if variable == 'pr':
        variable = 'pr_'

    """Find the paths to the directories containing the model data"""
    directory_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in directories:
            path = os.path.join(root, i)
            for j in list_of_models:
                if j in path and model_type in path and ensemble in path:
                    directory_paths = np.append(directory_paths, path)

    """Find the model files and their absolute paths."""
    model_file_paths = []
    for i in directory_paths:
        files = os.listdir(i)
        for j in files:
            if variable in j:
                model_file_path = os.path.join(i, j)
                model_file_paths = np.append(model_file_paths, model_file_path)

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    #model_file_paths = ['/ouce-home/data_not_backed_up/model/cmip5/CNRM-CM5/amip/mon/Lmon/r1i1p1/atlas_mrsos_Lmon_CNRM-CM5_amip_r1i1p1_197901-200812.nc']

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    """Load the data from the model file paths into a cube. Constrain the input years"""
    """Print the model ID, length of time dimension, and first and last model dates."""

    cubes = []

    if variable == 'hfls':
        name = 'surface_upward_latent_heat_flux'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'hfss':
        name = 'surface_upward_sensible_heat_flux'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'pr':
        name = 'precipitation_flux'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'mrsos':
        name = 'moisture_content_of_soil_layer'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    """Define a model number to begin with."""
    model_number = 0

    """For each cube (for each model),"""
    for model_data in cubes:

        """Select the model ID"""
        model_id = cubes[model_number].attributes['model_id']

        """If the variable is precipitation, multiply model cubes by 86400"""
        if variable == 'pr':
            model_data = iris.analysis.maths.multiply(cubes[model_number], 86400)

        """Reassign a model ID to the new multiplied cube."""
        model_data.long_name = model_id

        """Plot up a map for the file."""

        """ If the input month is defined as the whole year,"""
        if season_name == 'Climatology':

            """Take the mean over the cube."""
            model_data = model_data.collapsed('time', iris.analysis.MEAN)

        """ If the input month is defined as a season,"""
        if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

            """Create two coordinates on the cube to represent seasons and the season year."""
            iris.coord_categorisation.add_season(model_data, 'time', name='clim_season')
            iris.coord_categorisation.add_season_year(model_data, 'time', name='season_year')

            """Aggregate the data by season and season year."""
            seasonal_means = model_data.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

            """Constrain the data by the input season. Decapitalise the input variable."""
            constraint = iris.Constraint(clim_season=season_name.lower())
            model_data_season = seasonal_means.extract(constraint)

            """Take the mean over the cube."""
            model_data = model_data_season.collapsed('time', iris.analysis.MEAN)

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
        if variable == 'hfls':
            contour_levels = np.arange(80, 145, 5)
            cmap = matplotlib.cm.get_cmap('YlGnBu')
        if variable == 'hfss':
            contour_levels = np.arange(0, 65, 5)
            cmap = matplotlib.cm.get_cmap('YlGnBu_r')
        if variable == 'pr':
            contour_levels = np.arange(1, 12, 1)
            cmap = matplotlib.cm.get_cmap('YlGnBu')
        if variable == 'mrsos':
            contour_levels = np.arange(0, 65, 5)
            cmap = matplotlib.cm.get_cmap('YlGnBu')

        """Define the colour map and the projection."""
        crs_latlon = ccrs.PlateCarree()

        """Plot the map using cartopy, and add map features."""
        ax = plt.subplot(111, projection=crs_latlon)
        ax.set_extent([-22, 62, -22, 12], crs=crs_latlon)
        contour_plot = iplt.contourf(model_data, contour_levels, cmap=cmap, extend='both')
        ax.add_feature(lake_borders, zorder=5, edgecolor='k', linewidth=1)
        for i in country_borders:
            ax.add_geometries(i.geometry, ccrs.PlateCarree(), edgecolor="black", facecolor="None")
        ax.add_feature(cart.feature.OCEAN, zorder=1, facecolor="w")
        ax.add_feature(coastline, zorder=5, edgecolor='k', linewidth=2)

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

        if variable == 'hfls':
            colour_bar.set_ticks(np.arange(80, 145, 5))
            colour_bar.set_ticklabels(np.arange(80, 145, 5))
        if variable == 'hfss':
            colour_bar.set_ticks(np.arange(0, 65, 10))
            colour_bar.set_ticklabels(np.arange(0, 65, 10))
        if variable == 'pr':
            colour_bar.set_ticks(np.arange(1, 12, 1))
            colour_bar.set_ticklabels(np.arange(1, 12, 1))
        if variable == 'mrsos':
            colour_bar.set_ticks(np.arange(0, 65, 10))
            colour_bar.set_ticklabels(np.arange(0, 65, 10))

        colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)

        if variable_units == "W m-2":
            variable_units = "W $\mathregular{m^{-2}}$"
        if variable_units == 'kg m-2 s-1':
            variable_units = "mm $\mathregular{day^{-1}}$"
        if variable == 'mrsos':
            variable_name = "Volumetric Soil Moisture Content of Upper Layer"
            variable_units = "%"
        colour_bar.set_label(variable_name+" ("+variable_units+")", fontsize=10)

        """Add a title."""
        plt.title(model_data.long_name, fontsize=20)

        """Save the figure, close the plot and print an end statement."""
        fig.savefig(variable+"_"+season_name+"_"+model_id+".png")
        plt.close()
        print model_id+" plot done"

        """Add 1 to the model number to loop through the next model."""
        model_number +=1


#map(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], "amip", "hfss", "SON")

map(["IPSL-CM5B-LR/"], "amip", "mrsos", "SON")
