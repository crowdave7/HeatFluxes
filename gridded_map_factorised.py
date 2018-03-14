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
    start_time = time.time()
    model_file_paths = find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory)
    print time.time() - start_time, "seconds"

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
    print "cube_year_slicing"
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

def model_file_paths_ensemble(list_of_models, model_type, variable):
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


def contour_lev_colour_map(variable, lower_value, higher_value, interval, cmap):

    contour_levels = np.arange(lower_value, higher_value, interval)
    cmap = matplotlib.cm.get_cmap(cmap)
    return contour_levels, cmap


def colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, interval):

    colourbar_axis = fig.add_axes([0.20, 0.07, 0.60, 0.02])
    colour_bar = plt.colorbar(contour_plot, colourbar_axis, orientation='horizontal')
    colour_bar.set_ticks(np.arange(lower_tick, upper_tick, interval))
    colour_bar.set_ticklabels(np.arange(lower_tick, upper_tick, interval))
    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)

    return colour_bar

def plot_map(cubelist, variable, season_name, contour_levels, cmap, lower_tick, upper_tick, tick_interval):

    """Define variable name."""
    original_long_name = cubelist[0].long_name
    original_units = cubelist[0].units

    crs_latlon = ccrs.PlateCarree()

    fig = plt.figure(figsize=(8, 8))

    cube_number = 0

    for i in cubelist:
        long_name = i.long_name
        ax = plt.subplot(6, 4, cube_number+1, projection=crs_latlon)
        ax.set_extent([-22, 62, -24, 17], crs=crs_latlon)

        contour_plot = iplt.contourf(i, contour_levels, cmap=cmap, extend='both')

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

        if long_name == 'ACCESS1.3':
            long_name = 'ACCESS1-3'

        """Add a title."""
        plt.title(long_name, fontsize=8)

        """Add 1 to the model number to loop through the next model."""
        cube_number +=1
        print cube_number

    """Add colour bar."""
    colourbar_axis = fig.add_axes([0.20, 0.07, 0.60, 0.02])
    colour_bar = plt.colorbar(contour_plot, colourbar_axis, orientation='horizontal')

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

    fig.subplots_adjust(left=0.125, right=0.9, bottom=0.14, top=0.93, wspace=0.4, hspace=0.9)

    """Save the figure, close the plot and print an end statement."""
    print "saving final figure"
    fig.savefig("Gridded_"+variable+"_"+season_name+".png", dpi=600)
    plt.close()
    print "plot done"


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    model_type = "amip"
    list_of_reanalysis = ["cfsr", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    variable = "evapotranspiration"
    season_name = "SON"
    lower_year = 1979
    upper_year = 2008
    lower_value = 2.5
    higher_value = 5.25
    value_interval = 0.25
    lower_tick = 2.5
    upper_tick = 5.25
    tick_interval = 0.5
    cmap = "YlGnBu"
    time_range = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # MODELS AND ENSEMBLE MEAN

    if len(list_of_models) == 0:
        pass
    else:

        """Extract the model file paths."""
        model_file_paths = model_file_paths(list_of_models, model_type, variable)

        start_time = time.time()

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

        """Extract the regridded model file paths for the ensemble mean."""
        model_file_paths_ensemble = model_file_paths_ensemble(list_of_models, model_type, variable)

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
        ensemble_mean_cube.long_name = "Ensemble"

        """Append ensemble mean to the cubelist."""
        cubelist = np.append(cubelist, ensemble_mean_cube)

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # REANALYSIS

    if len(list_of_reanalysis) == 0:
        pass

    else:

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

        """Append reanalysis cubes to the cubelist."""
        cubelist = np.append(cubelist, cubelist_reanalysis)

    print cubelist

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # PLOTTING

    """Define the contour levels and colour map for the plot."""
    contour_levels = contour_lev_colour_map(variable, lower_value, higher_value, value_interval, cmap)[0]
    cmap = contour_lev_colour_map(variable, lower_value, higher_value, value_interval, cmap)[1]

    """Produce the plot."""
    plot_map(cubelist, variable, season_name, contour_levels, cmap, lower_tick, upper_tick, tick_interval)
