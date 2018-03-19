#!/usr/bin/env python

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
import multiprocessing
import numpy as np
import os
import time



def build_cubelist(i, array):

    cubes_in_file = iris.load(model_file_paths[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        if len(coord_names) == 3:
            array[i] = cube

def model_file_paths(list_of_models, model_type, variable):

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

    """FIND MODEL PATHS."""
    print "finding model paths"
    model_file_paths = find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory)

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """If variable is evspsbl_, convert variable back to evspsbl"""
    if variable == 'evspsbl_':
        variable = 'evspsbl'

    return model_file_paths

def find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory):
    "model_file_paths"

    if variable == 'evapotranspiration':
        variable = 'hfls'

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
    return model_file_paths

def run_cube_year_slicing(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        if variable != "treeFrac":
            cubelist[i] = cubelist[i].extract(time_range)
            time_points = cubelist[i].coord('time').points
            times = cubelist[i].coord('time').units.num2date(time_points)
        model_id = cubelist[i].attributes['model_id']
        print model_id
        if variable != "treeFrac":
            print len(times)
            print times[0]
            print times[-1]

    array[i] = cubelist[i]

def run_cube_season_slicing_change_units(i, array):

    """Slice the regridded cube down to the African domain."""
    cubelist[i] = cubelist[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

    """Select model ID."""
    if variable == 'nrad':
        model_id = cubelist[i].long_name
    else:
        model_id = cubelist[i].attributes['model_id']

    if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg':
        cubelist[i] = iris.analysis.maths.multiply(cubelist[i], 86400)

    """if variable is evapotranspiration, divide by 28"""
    if variable == 'evapotranspiration':
        cubelist[i] = iris.analysis.maths.divide(cubelist[i], 28)

    """Reassign model ID."""
    cubelist[i].long_name = model_id

    """ If the input month is defined as the whole year,"""
    if season_name == 'Climatology':

        """Take the mean over the cube."""
        cubelist[i] = cubelist[i].collapsed('time', iris.analysis.MEAN)

    """ If the input month is defined as a season,"""
    if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

        if variable == 'treeFrac':
            cubelist[i] = cubelist[i].collapsed('time', iris.analysis.MEAN)

        else:
            """Create two coordinates on the cube to represent seasons and the season year."""
            iris.coord_categorisation.add_season(cubelist[i], 'time', name='clim_season')
            iris.coord_categorisation.add_season_year(cubelist[i], 'time', name='season_year')
            """Aggregate the data by season and season year."""
            seasonal_means = cubelist[i].aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)
            """Constrain the data by the input season. Decapitalise the input variable."""
            constraint = iris.Constraint(clim_season=season_name.lower())
            cube_season = seasonal_means.extract(constraint)
            """Take the mean over the cube."""
            cubelist[i] = cube_season.collapsed('time', iris.analysis.MEAN)

    array[i] = cubelist[i]

def build_cubelist_ensemble(i, array):

    cubes_in_file = iris.load(model_file_paths_ensemble[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        if len(coord_names) == 3:
            array[i] = cube

def model_file_paths_ensemble_func(list_of_models, model_type, variable):

        """Import the data."""
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles"

        """If variable is pr, distinguish between pr and precipitable water to find model files."""
        if variable == 'pr':
            variable = 'pr_'

        """If variable is evspsbl, distinguish between evspsbl and evspsblsoi to find model files."""
        if variable == 'evspsbl':
            variable = 'evspsbl_'

        if variable == 'evapotranspiration':
            variable = 'hfls'

        model_file_paths = find_model_paths_ensemble(list_of_models, model_type, variable, root_directory)

        """If variable is pr_, convert variable back to pr"""
        if variable == 'pr_':
            variable = 'pr'

        """If variable is evspsbl_, convert variable back to evspsbl"""
        if variable == 'evspsbl_':
            variable = 'evspsbl'

        if variable == 'hfls':
            variable = 'evapotranspiration'

        return model_file_paths

def find_model_paths_ensemble(list_of_models, model_type, variable, root_directory):

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

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    return model_file_paths

def run_cube_year_slicing_ensemble(i, array):
    print "cube_year_slicing"
    with iris.FUTURE.context(cell_datetime_objects=True):
        if variable != "treeFrac":
            cubelist_ensemble[i] = cubelist_ensemble[i].extract(time_range)
            time_points = cubelist_ensemble[i].coord('time').points
            times = cubelist_ensemble[i].coord('time').units.num2date(time_points)
        model_id = cubelist_ensemble[i].attributes['model_id']
        print model_id
        if variable != "treeFrac":
            print len(times)
            print times[0]
            print times[-1]

    array[i] = cubelist_ensemble[i]

def run_cube_season_slicing_change_units_ensemble(i, array):

    """Slice the regridded cube down to the African domain."""
    cubelist_ensemble[i] = cubelist_ensemble[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

    """Select model ID."""
    if variable == 'nrad':
        model_id = cubelist_ensemble[i].long_name
    else:
        model_id = cubelist_ensemble[i].attributes['model_id']

    if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg':
        cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 86400)

    """if variable is evapotranspiration, divide by 28"""
    if variable == 'evapotranspiration':
        cubelist_ensemble[i] = iris.analysis.maths.divide(cubelist_ensemble[i], 28)

    """Reassign model ID."""
    cubelist_ensemble[i].long_name = model_id

    """ If the input month is defined as the whole year,"""
    if season_name == 'Climatology':

        """Take the mean over the cube."""
        cubelist_ensemble[i] = cubelist_ensemble[i].collapsed('time', iris.analysis.MEAN)

    """ If the input month is defined as a season,"""
    if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

        if variable == 'treeFrac':
            cubelist_ensemble[i] = cubelist_ensemble[i].collapsed('time', iris.analysis.MEAN)

        else:
            """Create two coordinates on the cube to represent seasons and the season year."""
            iris.coord_categorisation.add_season(cubelist_ensemble[i], 'time', name='clim_season')
            iris.coord_categorisation.add_season_year(cubelist_ensemble[i], 'time', name='season_year')

            """Aggregate the data by season and season year."""
            seasonal_means = cubelist_ensemble[i].aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

            """Constrain the data by the input season. Decapitalise the input variable."""
            constraint = iris.Constraint(clim_season=season_name.lower())
            cube_season = seasonal_means.extract(constraint)

            """Take the mean over the cube."""
            cubelist_ensemble[i] = cube_season.collapsed('time', iris.analysis.MEAN)

    array[i] = cubelist_ensemble[i]

def build_cubelist_reanalysis(i, array):

    cubes_in_file = iris.load(reanalysis_file_paths[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        if len(coord_names) == 3:
            array[i] = cube

def reanalysis_file_paths(list_of_reanalysis, variable):
    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles"

    """If variable is pr, distinguish between pr and precipitable water to find model files."""
    if variable == 'pr':
        variable = 'pr_'

    """If variable is evspsbl, distinguish between evspsbl and evspsblsoi to find model files."""
    if variable == 'evspsbl':
        variable = 'evspsbl_'

    if variable == 'evapotranspiration':
        variable = 'hfls'

    reanalysis_file_paths = find_reanalysis_file_paths(list_of_reanalysis, variable, root_directory)

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """If variable is evspsbl_, convert variable back to evspsbl"""
    if variable == 'evspsbl_':
        variable = 'evspsbl'

    if variable == 'hfls':
        variable = 'evapotranspiration'

    return reanalysis_file_paths

def find_reanalysis_file_paths(list_of_reanalysis, variable, root_directory):

    """Find the paths to the files containing the reanalysis data"""
    reanalysis_file_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            for j in list_of_reanalysis:
                if j in path and variable in path:
                    reanalysis_file_paths = np.append(reanalysis_file_paths, path)

    reanalysis_file_paths = sorted(reanalysis_file_paths, key=lambda s: s.lower())

    return reanalysis_file_paths


def run_cube_year_slicing_reanalysis(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        if variable != "treeFrac":
            cubelist_reanalysis[i] = cubelist_reanalysis[i].extract(time_range)
            time_points = cubelist_reanalysis[i].coord('time').points
            times = cubelist_reanalysis[i].coord('time').units.num2date(time_points)
        reanalysis_id = list_of_reanalysis[i]
        print reanalysis_id
        if variable != "treeFrac":
            print len(times)
            print times[0]
            print times[-1]

    array[i] = cubelist_reanalysis[i]


def run_cube_season_slicing_change_units_reanalysis(i, array):

    """Slice the regridded cube down to the African domain."""
    cubelist_reanalysis[i] = cubelist_reanalysis[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

    """if variable is evapotranspiration, divide by 28"""
    if variable == 'evapotranspiration':
        cubelist_reanalysis[i] = iris.analysis.maths.divide(cubelist_reanalysis[i], 28)

    """ If the input month is defined as the whole year,"""
    if season_name == 'Climatology':

        """Take the mean over the cube."""
        cubelist_reanalysis[i] = cubelist_reanalysis[i].collapsed('time', iris.analysis.MEAN)

    """ If the input month is defined as a season,"""
    if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

        """Create two coordinates on the cube to represent seasons and the season year."""
        iris.coord_categorisation.add_season(cubelist_reanalysis[i], 'time', name='clim_season')
        iris.coord_categorisation.add_season_year(cubelist_reanalysis[i], 'time', name='season_year')

        """Aggregate the data by season and season year."""
        seasonal_means = cubelist_reanalysis[i].aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

        """Constrain the data by the input season. Decapitalise the input variable."""
        constraint = iris.Constraint(clim_season=season_name.lower())
        reanalysis_data_season = seasonal_means.extract(constraint)

        """Take the mean over the cube."""
        cubelist_reanalysis[i] = reanalysis_data_season.collapsed('time', iris.analysis.MEAN)

        """Select the reanalysis ID."""
        reanalysis_id = list_of_reanalysis[i]

        if reanalysis_id == "cfsr":
            cubelist_reanalysis[i].long_name = "CFSR"
            cubelist_reanalysis[i].rename("CFSR")
        if reanalysis_id == "erai":
            cubelist_reanalysis[i].long_name = "ERA-Interim"
            cubelist_reanalysis[i].rename("ERA-Interim")
        if reanalysis_id == "gleam":
            cubelist_reanalysis[i].long_name = "GLEAM-LE"
            cubelist_reanalysis[i].rename("GLEAM-LE")
        if reanalysis_id == "jra":
            cubelist_reanalysis[i].long_name = "JRA-55"
            cubelist_reanalysis[i].rename("JRA-55")
        if reanalysis_id == "merra2":
            cubelist_reanalysis[i].long_name = "MERRA-2"
            cubelist_reanalysis[i].rename("MERRA-2")
        if reanalysis_id == "mswep":
            cubelist_reanalysis[i].long_name = "GLEAM (MSWEP)"
            cubelist_reanalysis[i].rename("GLEAM (MSWEP)")
        if reanalysis_id == "ncep-ncar":
            cubelist_reanalysis[i].long_name = "NCEP/NCAR"
            cubelist_reanalysis[i].rename("NCEP/NCAR")
        if reanalysis_id == "ncep-doe":
            cubelist_reanalysis[i].long_name = "NCEP DOE-2"
            cubelist_reanalysis[i].rename("NCEP DOE-2")

    array[i] = cubelist_reanalysis[i]

def contour_lev_colour_map(lower_value, higher_value, interval, cmap):

    contour_levels = np.arange(lower_value, higher_value, interval)
    cmap = matplotlib.cm.get_cmap(cmap)
    return contour_levels, cmap

def colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, interval):

    colour_bar.set_ticks(np.arange(lower_tick, upper_tick, interval))
    colour_bar.set_ticklabels(np.arange(lower_tick, upper_tick, interval))
    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)

    return colour_bar

def plot_map(cube, variable, season_name, contour_levels, cmap, lower_tick, upper_tick, tick_interval):

    """Define variable name."""
    original_long_name = cube.long_name
    original_units = cube.units

    crs_latlon = ccrs.PlateCarree()

    """Plot the figure."""
    fig = plt.figure()

    ax = plt.subplot(111, projection=crs_latlon)
    #ax.set_extent([-22, 62, -22, 12], crs=crs_latlon)

    contour_plot = iplt.contourf(cube, contour_levels, cmap=cmap, extend='both')

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
    #ax.add_feature(cart.feature.OCEAN, zorder=1, facecolor="w")

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

    if cube.long_name == 'ACCESS1.3':
        cube.long_name = 'ACCESS1-3'

    """Add a title."""
    plt.title(cube.long_name, fontsize=14)

    """Add a colour bar, with ticks and labels."""
    colour_bar = plt.colorbar(contour_plot, orientation='horizontal', pad=0.1, aspect=40)

    """Adjust ticks."""
    colour_bar = colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, tick_interval)

    """Redefine variable name."""
    variable_name = str(original_long_name)
    variable_units = str(original_units)

    if variable_units == "W m-2":
        variable_units = "W $\mathregular{m^{-2}}$"
    if variable_units == 'kg m-2 s-1':
        variable_units = "mm $\mathregular{day^{-1}}$"
    if variable == 'mrsos':
        variable_name = "Soil Moisture Content of Upper Layer"
        variable_units = "mm"
    if variable == 'nrad':
        variable_name = "Surface Net Radiation"
        variable_units = "W $\mathregular{m^{-2}}$"
    if variable == 'evapotranspiration':
        variable_name = 'Evapotranspiration'
        variable_units = "mm $\mathregular{day^{-1}}$"
    if variable == 'nrad':
        variable_name = "Surface Net Radiation"
        variable_units = "W $\mathregular{m^{-2}}$"

    colour_bar.set_label(variable_name+" ("+variable_units+")", fontsize=10)

    cube_name = cube.long_name

    if cube_name == "Composite":
        """Save the figure, close the plot and print an end statement."""
        print "hi1"
        fig.savefig(variable+"_"+season_name+"_Composite.png")
        plt.close()
        print "composite plot done"

    if cube_name == "Composite Difference":
        """Save the figure, close the plot and print an end statement."""
        fig.savefig(variable+"_"+season_name+"_Composite_Difference.png")
        plt.close()
        print "composite difference plot done"

    else:
        print "hi3"
        fig.savefig(variable+"_"+season_name+"_"+cube.long_name+".png")
        plt.close()
        print cube.long_name+" plot done"


if __name__ == "__main__":

    # CANNOT PRINT COMPOSITE MAPS OF BOTH MODEL AND REANALYSIS AT ONCE.
    # CANNOT HAVE MORE THAN 1 REANALYSIS.
    list_of_models = ['inmcm4/', 'MIROC5/', 'BNU-ESM/', 'ACCESS1-3/', 'GFDL-HIRAM-C360/']
    #list_of_models = ['MRI-AGCM3-2H/', 'bcc-csm1-1/', 'CanAM4/', 'bcc-csm1-1-m/', 'MRI-CGCM3/']
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    list_of_reanalysis = []

    list_of_models_higher_composite = ['inmcm4/', 'MIROC5/', 'BNU-ESM/', 'ACCESS1-3/', 'GFDL-HIRAM-C360/']
    list_of_models_lower_composite = ['MRI-AGCM3-2H/', 'bcc-csm1-1/', 'CanAM4/', 'bcc-csm1-1-m/', 'MRI-CGCM3/']

    #list_of_reanalysis = ["cfsr", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    model_type = "amip"
    variable = "evapotranspiration"
    season_name = "SON"
    lower_year = 1979
    upper_year = 2008

    lower_value = 2.5
    higher_value = 5.75
    value_interval = 0.25
    lower_tick = 2.5
    upper_tick = 5.75
    tick_interval = 0.5

    lower_value_diff = -2.0
    higher_value_diff = 2.25
    value_interval_diff = 0.25
    lower_tick_diff = -2.0
    upper_tick_diff = 2.25
    tick_interval_diff = 0.5

    composite_difference = "yes"

    if composite_difference == "yes":
        cmap = "bwr"
    else:
        cmap = "YlGnBu"

    time_range = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)


    if composite_difference != "yes":

        # ------------------------------------------------------------------------------------------------------------------------------------------
        # MODELS AND ENSEMBLE MEAN

        if len(list_of_models) == 0:
            pass

        if len(list_of_models) == 1:
            """Extract the model file paths."""
            start_time = time.time()
            model_file_paths = model_file_paths(list_of_models, model_type, variable)
            print time.time() - start_time, "seconds"

            """Build a list of cubes from the model file paths."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(model_file_paths)):
                p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            """Define original long name and original units for later."""
            original_long_name = cubelist[0].long_name
            original_unit = cubelist[0].units

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist)):
                p = multiprocessing.Process(target=run_cube_year_slicing, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            """Slice each cube by the season name, and change units if required."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist)):
                p = multiprocessing.Process(target=run_cube_season_slicing_change_units, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            model_cube = cubelist[0]

            contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
            cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

            plot_map(model_cube, variable, season_name, contour_levels, cmap, lower_tick, upper_tick, tick_interval)

        if len(list_of_models) > 1:

            """Extract the regridded model file paths for the ensemble mean."""
            print list_of_models
            model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

            """Build a list of cubes from the regridded model file paths."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(model_file_paths_ensemble)):
                p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=run_cube_year_slicing_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Slice each cube by the season name, and change units if required."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=run_cube_season_slicing_change_units_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Compute the ensemble mean cube."""
            ensemble_mean_cube = sum(cubelist_ensemble) / float(len(cubelist_ensemble))
            ensemble_mean_cube.long_name = "Composite"

            contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
            cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

            plot_map(ensemble_mean_cube, variable, season_name, contour_levels, cmap, lower_tick, upper_tick, tick_interval)

        if len(list_of_reanalysis) == 0:
            pass

        if len(list_of_reanalysis) == 1:

            """Build a list of cubes from the reanalysis file paths."""
            reanalysis_file_paths = reanalysis_file_paths(list_of_reanalysis, variable)
            """Load the cubes."""

            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(reanalysis_file_paths)):
                p = multiprocessing.Process(target=build_cubelist_reanalysis, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_reanalysis = array.values()

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_reanalysis)):
                p = multiprocessing.Process(target=run_cube_year_slicing_reanalysis, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_reanalysis = array.values()

            """Slice each cube by the season name, and change units if required."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_reanalysis)):
                p = multiprocessing.Process(target=run_cube_season_slicing_change_units_reanalysis, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_reanalysis = array.values()

            reanalysis_cube = cubelist_reanalysis[0]

            contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
            cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

            plot_map(reanalysis_cube, variable, season_name, contour_levels, cmap, lower_tick, upper_tick, tick_interval)

    if composite_difference == "yes":

        if len(list_of_models_higher_composite) == 0:
            pass

        if len(list_of_models_higher_composite) == 1:

            """Extract the model file paths."""
            start_time = time.time()
            model_file_paths = model_file_paths(list_of_models_higher_composite, model_type, variable)
            print time.time() - start_time, "seconds"

            """Build a list of cubes from the model file paths."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(model_file_paths)):
                p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            """Define original long name and original units for later."""
            original_long_name = cubelist[0].long_name
            original_unit = cubelist[0].units

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist)):
                p = multiprocessing.Process(target=run_cube_year_slicing, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            """Slice each cube by the season name, and change units if required."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist)):
                p = multiprocessing.Process(target=run_cube_season_slicing_change_units, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            higher_composite_cube = cubelist[0]

        if len(list_of_models_higher_composite) > 1:

            """Extract the regridded model file paths for the ensemble mean."""
            model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models_higher_composite, model_type, variable)

            print model_file_paths_ensemble

            """Build a list of cubes from the regridded model file paths."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(model_file_paths_ensemble)):
                p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=run_cube_year_slicing_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Slice each cube by the season name, and change units if required."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=run_cube_season_slicing_change_units_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Compute the ensemble mean cube."""
            higher_composite_cube = sum(cubelist_ensemble) / float(len(cubelist_ensemble))
            higher_composite_cube.long_name = "Strong Composite"

            print higher_composite_cube

        if len(list_of_models_lower_composite) == 0:
            pass

        if len(list_of_models_lower_composite) == 1:
            """Extract the model file paths."""

            start_time = time.time()
            model_file_paths = model_file_paths(list_of_models_lower_composite, model_type, variable)
            print time.time() - start_time, "seconds"

            """Build a list of cubes from the model file paths."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(model_file_paths)):
                p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            """Define original long name and original units for later."""
            original_long_name = cubelist[0].long_name
            original_unit = cubelist[0].units

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist)):
                p = multiprocessing.Process(target=run_cube_year_slicing, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            """Slice each cube by the season name, and change units if required."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist)):
                p = multiprocessing.Process(target=run_cube_season_slicing_change_units, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist = array.values()

            lower_composite_cube = cubelist[0]

        if len(list_of_models_lower_composite) > 1:

            print list_of_models
            print list_of_models_lower_composite

            """Extract the regridded model file paths for the ensemble mean."""
            model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models_lower_composite, model_type, variable)

            """Build a list of cubes from the regridded model file paths."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(model_file_paths_ensemble)):
                p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=run_cube_year_slicing_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Slice each cube by the season name, and change units if required."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=run_cube_season_slicing_change_units_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Compute the ensemble mean cube."""
            lower_composite_cube = sum(cubelist_ensemble) / float(len(cubelist_ensemble))
            lower_composite_cube.long_name = "Weak Composite"

            print lower_composite_cube

        composite_difference_cube = higher_composite_cube - lower_composite_cube
        composite_difference_cube.long_name = "Composite Difference"

        contour_levels = contour_lev_colour_map(lower_value_diff, higher_value_diff, value_interval_diff, cmap)[0]
        cmap = contour_lev_colour_map(lower_value_diff, higher_value_diff, value_interval_diff, cmap)[1]

        plot_map(composite_difference_cube, variable, season_name, contour_levels, cmap, lower_tick_diff, upper_tick_diff, tick_interval_diff)
