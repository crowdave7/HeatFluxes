#!/usr/bin/env python

"""Import necessary modules for this code."""
import cartopy as cart
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import ensemble_cube
import iris
import iris.analysis
import iris.coord_categorisation
import iris.plot as iplt
import matplotlib
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np
import os
import reanalysis_cube

# WHEN ADDING NEW VARIABLE, MODIFY THIS SCRIPT, REANALYSIS_CUBE AND ENSEMBLE_CUBE AND REGRID NEW VARIABLES

def map(list_of_models, model_type, variable, season_name):
    """Take the input variables, and find the paths to the relevant model files."""
    """Plot these up as a spatial map for the given season."""
    """Central African domain."""

    """Import the data."""

    if variable == 'nrad':
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1NetRadiationModelFiles"
    else:
        root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

    """If variable is pr, distinguish between pr and precipitable water to find model files."""
    if variable == 'pr':
        variable = 'pr_'

    """If variable is evspsbl, distinguish between evspsbl and evspsblsoi to find model files."""
    if variable == 'evspsbl':
        variable = 'evspsbl_'

    variable_original = []

    if variable == 'evapotranspiration':
        variable = 'hfls'
        variable_original = 'evapotranspiration'

    if variable == 'nrad':
        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                for j in list_of_models:
                    if j == "bcc-csm1-1/":
                        j = "bcc-csm1-1_"
                    for char in '/':
                        j = j.replace(char,'')
                    if j in path and model_type in path and variable in path:
                        model_file_paths = np.append(model_file_paths, path)
    else:

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
    print model_file_paths

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """If variable is evspsbl_, convert variable back to evspsbl"""
    if variable == 'evspsbl_':
        variable = 'evspsbl'

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

    if variable == 'nrad':
        for i in model_file_paths:
            cube = iris.load_cube(i)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = i.long_name
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print len(times)
                print model_id
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

    if variable == 'tran':
        name = 'transpiration_flux'
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

    if variable == 'evspsblsoi':
        name = 'water_evaporation_flux_from_soil'
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

    if variable == 'evspsbl':
        name = 'water_evaporation_flux'
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

    if variable == 'evspsblveg':
        name = 'water_evaporation_flux_from_canopy'
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

    if variable == 'treeFrac':
        for i in model_file_paths:
            cube = iris.load_cube(i)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                count +=1

    if variable == 'lai':
        for i in model_file_paths:
            name = 'leaf_area_index'
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

    """Define the contour levels for the input variables."""
    if variable == 'hfls' and variable_original != 'evapotranspiration':
        contour_levels = np.arange(80, 145, 5)
    if variable == 'hfls' and variable_original == 'evapotranspiration':
        contour_levels = np.arange(2.5, 5.25, 0.25)
    if variable == 'hfss':
        contour_levels = np.arange(0, 65, 5)
    if variable == 'pr':
        contour_levels = np.arange(1, 12, 1)
    if variable == 'mrsos':
        contour_levels = np.arange(10, 55, 3)
    if variable == 'tran':
        contour_levels = np.arange(0.2, 3.2, 0.2)
    if variable == 'evspsblsoi':
        contour_levels = np.arange(0.2, 2.2, 0.2)
    if variable == 'evspsbl':
        contour_levels = np.arange(2.0, 5.2, 0.2)
    if variable == 'evspsblveg':
        contour_levels = np.arange(0.2, 2.2, 0.2)
    if variable == 'nrad':
        contour_levels = np.arange(100, 185, 5)
    if variable == 'treeFrac':
        contour_levels = np.arange(0, 110, 10)
    # if variable == 'lai':
    #     contour_levels = np.arange(0, 20, 1)

    """Define the colour map and the projection."""
    if variable == 'hfls' and variable_original != 'evapotranspiration':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'hfls' and variable_original == 'evapotranspiration':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'hfss':
        cmap = matplotlib.cm.get_cmap('YlGnBu_r')
    if variable == 'pr':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'nrad':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'mrsos':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'tran':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'evspsblsoi':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'evspsbl':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'evspsblveg':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'treeFrac':
        cmap = matplotlib.cm.get_cmap('YlGnBu')
    if variable == 'lai':
        cmap = matplotlib.cm.get_cmap('YlGnBu')

    crs_latlon = ccrs.PlateCarree()

    """Define a model number to begin with."""
    model_number = 0

    """Define the original long name."""
    original_long_name = cubes[0].long_name
    original_units = cubes[0].units

    """Plot the figure."""
    fig = plt.figure(figsize=(8, 8))

    """For each cube (for each model),"""
    for model_data in cubes:

        """Select the model ID"""
        if variable == 'nrad':
            model_id = model_data.long_name
        else:
            model_id = cubes[model_number].attributes['model_id']

        """If the variable is precipitation, or land surface quantity (kg m-2 s-1) multiply model cubes by 86400"""
        if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg':
            model_data = iris.analysis.maths.multiply(cubes[model_number], 86400)

        """if variable is evapotranspiration, divide by 28"""
        if variable_original == 'evapotranspiration':
            model_data = iris.analysis.maths.divide(cubes[model_number], 28)

        """Reassign a model ID to the new multiplied cube."""
        model_data.long_name = model_id

        """Plot up a map for the file."""

        """ If the input month is defined as the whole year,"""
        if season_name == 'Climatology':

            """Take the mean over the cube."""
            model_data = model_data.collapsed('time', iris.analysis.MEAN)

            """Add 1 to the model number to loop through the next model."""
            model_number +=1

        """ If the input month is defined as a season,"""
        if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

            if variable == 'treeFrac':
                model_data = model_data.collapsed('time', iris.analysis.MEAN)

            else:

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

            """Return this mean to the cubelist."""
            cubes[model_number] = model_data

            print model_number
            """Add 1 to the model number to loop through the next model."""
            model_number +=1

    if variable_original == 'evapotranspiration':
        variable = 'evapotranspiration'

    """Add the ensemble mean to the cubelist."""
    ensemble_cube_data = ensemble_cube.ensemble(list_of_models, model_type, variable, season_name)
    cubes = np.append(cubes, ensemble_cube_data)

    """Ignore reanalysis if dealing with partitioned evapotranspiration variables."""

    if variable == 'evspsbl' or variable == 'evspsblsoi' or variable == 'evspsblveg' or variable == 'tran' or variable == 'nrad' or variable == 'treeFrac' or variable == 'lai':
        pass

    else:

        if variable != 'tran' or variable != 'evspsblsoi' or variable != 'evspsbl' or variable != 'evspsblveg':
            """Add CFSR data to the cubelist."""
            cfsr_cube_data = reanalysis_cube.reanalysis(["cfsr"], variable, season_name)
            cubes = np.append(cubes, cfsr_cube_data)

        if variable != 'tran' or variable != 'evspsblsoi' or variable != 'evspsbl' or variable != 'evspsblveg':
            """Add ERA-Interim data to the cubelist."""
            erai_cube_data = reanalysis_cube.reanalysis(["erai"], variable, season_name)
            cubes = np.append(cubes, erai_cube_data)

        if variable == 'hfss':
            """Add GLEAM data to the cubelist."""
            gleam_cube_data = reanalysis_cube.reanalysis(["gleam"], variable, season_name)
            cubes = np.append(cubes, gleam_cube_data)

        if variable == 'hfls':
            """Add GLEAM data to the cubelist."""
            gleam_cube_data = reanalysis_cube.reanalysis(["gleam"], variable, season_name)
            cubes = np.append(cubes, gleam_cube_data)

        if variable == 'evapotranspiration':
            """Add GLEAM data to the cubelist."""
            gleam_cube_data = reanalysis_cube.reanalysis(["gleam"], variable, season_name)
            cubes = np.append(cubes, gleam_cube_data)

        if variable == 'mrsos':
            """Add GLEAM data to the cubelist."""
            gleam_cube_data = reanalysis_cube.reanalysis(["gleam"], variable, season_name)
            cubes = np.append(cubes, gleam_cube_data)

        if variable == 'pr':
            """Add MSWEP data to the cubelist."""
            mswep_cube_data = reanalysis_cube.reanalysis(["mswep"], variable, season_name)
            cubes = np.append(cubes, mswep_cube_data)

        if variable != "mrsos":
        #if variable != "mrsos" or variable != 'tran' or variable != 'evspsblsoi' or variable != 'evspsbl':
            """Add JRA-55 data to the cubelist."""
            jra_cube_data = reanalysis_cube.reanalysis(["jra"], variable, season_name)
            cubes = np.append(cubes, jra_cube_data)

        if variable != 'tran' or variable != 'evspsblsoi' or variable != 'evspsbl':
            """Add MERRA-2 data to the cubelist."""
            merra_cube_data = reanalysis_cube.reanalysis(["merra2"], variable, season_name)
            cubes = np.append(cubes, merra_cube_data)

        if variable != 'tran' or variable != 'evspsblsoi' or variable != 'evspsbl':
            """Add NCEP DOE-2 data to the cubelist."""
            doe_cube_data = reanalysis_cube.reanalysis(["ncep-doe"], variable, season_name)
            cubes = np.append(cubes, doe_cube_data)

    #print doe_cube_data
    print cubes

    """Plot each model and the ensemble mean up."""
    model_number = 0

    for model_data in cubes:

        ax = plt.subplot(6, 4, model_number+1, projection=crs_latlon)
        ax.set_extent([-22, 62, -24, 17], crs=crs_latlon)
        contour_plot = iplt.contourf(model_data, cmap=cmap, extend='both')

        """Import coastlines and lake borders. Set the scale to 10m, 50m or 110m resolution for more detail."""
        coastline = cart.feature.NaturalEarthFeature(category='physical', name='coastline', scale='110m', facecolor='none')
        lake_borders = cart.feature.NaturalEarthFeature(category='physical', name='lakes', scale='110m', facecolor='none')

        """Import country borders."""
        shapefile = shapereader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries')
        reader = shapereader.Reader(shapefile)
        country_borders = reader.records()

        """Remove iris warning message."""
        iris.FUTURE.netcdf_promote = True

        """Plot the map using cartopy, and add map features."""

        ax.add_feature(coastline, zorder=5, edgecolor='k', linewidth=2)
        ax.add_feature(lake_borders, zorder=5, edgecolor='k', linewidth=1)
        for i in country_borders:
            ax.add_geometries(i.geometry, ccrs.PlateCarree(), edgecolor="black", facecolor="None")
        ax.add_feature(cart.feature.OCEAN, zorder=1, facecolor="w")

        """Define gridlines."""
        gridlines = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='black', linewidth=0.4, linestyle='--')
        gridlines.xlabels_top = False
        gridlines.xlabels_bottom = True
        gridlines.ylabels_left = True
        gridlines.ylabels_right = False
        gridlines.xlabel_style = {'size': 6, 'color': 'black'}
        gridlines.ylabel_style = {'size': 6, 'color': 'black'}
        gridlines.xlines = False
        gridlines.ylines = False
        gridlines.xformatter = LONGITUDE_FORMATTER
        gridlines.yformatter = LATITUDE_FORMATTER
        gridlines.xlocator = mticker.FixedLocator(np.arange(-40, 100, 20))
        gridlines.ylocator = mticker.FixedLocator(np.arange(-50, 70, 10))

        if model_data.long_name == 'ACCESS1.3':
            model_data.long_name = 'ACCESS1-3'

        """Add a title."""
        plt.title(model_data.long_name, fontsize=8)

        """Add 1 to the model number to loop through the next model."""
        model_number +=1
        print model_number

    colourbar_axis = fig.add_axes([0.20, 0.07, 0.60, 0.02])
    colour_bar = plt.colorbar(contour_plot, colourbar_axis, orientation='horizontal')

    if variable == 'hfls' and variable_original != 'evapotranspiration':
        colour_bar.set_ticks(np.arange(80, 145, 10))
        colour_bar.set_ticklabels(np.arange(80, 145, 10))
    if variable == 'hfls' and variable_original == 'evapotranspiration':
        colour_bar.set_ticks(np.arange(2.5, 5.25, 0.5))
        colour_bar.set_ticklabels(np.arange(2.5, 5.25, 0.5))
    if variable == 'hfss':
        colour_bar.set_ticks(np.arange(0, 65, 10))
        colour_bar.set_ticklabels(np.arange(0, 65, 10))
    if variable == 'pr':
        colour_bar.set_ticks(np.arange(1, 12, 1))
        colour_bar.set_ticklabels(np.arange(1, 12, 1))
    if variable == 'nrad':
        colour_bar.set_ticks(np.arange(100, 185, 10))
        colour_bar.set_ticklabels(np.arange(100, 185, 10))
    if variable == 'mrsos':
        colour_bar.set_ticks(np.arange(10, 55, 6))
        colour_bar.set_ticklabels(np.arange(10, 55, 6))
    if variable == 'tran':
        colour_bar.set_ticks(np.arange(0.2, 3.2, 0.2))
        colour_bar.set_ticklabels(np.arange(0.2, 3.2, 0.2))
    if variable == 'evspsblsoi':
        colour_bar.set_ticks(np.arange(0.2, 2.2, 0.2))
        colour_bar.set_ticklabels(np.arange(0.2, 2.2, 0.2))
    if variable == 'evspsbl':
        colour_bar.set_ticks(np.arange(2.0, 5.2, 0.2))
        colour_bar.set_ticklabels(np.arange(2.0, 5.2, 0.2))
    if variable == 'evspsveg':
        colour_bar.set_ticks(np.arange(0.2, 2.2, 0.2))
        colour_bar.set_ticklabels(np.arange(0.2, 2.2, 0.2))
    if variable == 'evspsveg':
        colour_bar.set_ticks(np.arange(0, 110, 10))
        colour_bar.set_ticklabels(np.arange(0, 110, 10))

    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)

    variable_name = str(original_long_name)
    variable_units = str(original_units)

    if variable_units == "W m-2":
        variable_units = "W $\mathregular{m^{-2}}$"
    if variable_units == 'kg m-2 s-1':
        variable_units = "mm $\mathregular{day^{-1}}$"
    if variable == 'mrsos':
        variable_name = "Soil Moisture Content of Upper Layer"
        variable_units = "mm"
    if variable_original == 'evapotranspiration':
        variable_name = 'Evapotranspiration'
        variable_units = "mm $\mathregular{day^{-1}}$"
    if variable == 'nrad':
        variable_name = "Surface Net Radiation"
        variable_units = "W $\mathregular{m^{-2}}$"

    colour_bar.set_label(variable_name+" ("+variable_units+")", fontsize=10)

    fig.subplots_adjust(left=0.125, right=0.9, bottom=0.14, top=0.93, wspace=0.4, hspace=0.9)

    """Save the figure, close the plot and print an end statement."""
    print "saving final figure"
    fig.savefig("Gridded_"+variable+"_"+season_name+".png", dpi=600)
    plt.close()
    print "plot done"

#map(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], "amip", "mrsos", "SON")

#map(["ACCESS1-0/", "ACCESS1-3/", "BNU-ESM/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR"], "amip", "treeFrac", "SON")

map(["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1_", "bcc-csm1-1-m", "BNU-ESM/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "inmcm4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M"], "amip", "lai", "SON")



#map(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], "amip", "mrsos", "SON")

#map(["bcc-csm1-1/", "BNU-ESM", "CanAM4", "GFDL-HIRAM-C360", "GISS-E2-R/", "inmcm4", "MIROC5", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], "amip", "evspsbl", "SON")

#map(["bcc-csm1-1/", "BNU-ESM", "CanAM4", "GFDL-HIRAM-C360", "GISS-E2-R/", "inmcm4", "MIROC5", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], "amip", "tran", "SON")


#map(["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "nrad", "SON")

#map(["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "mrsos", "SON")

#map(["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "tran", "SON")
