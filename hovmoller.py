#!/usr/bin/env python

"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.analysis.cartography
import iris.coord_categorisation
import iris.util
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans, interp
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
import time
from netCDF4 import num2date
from netCDF4 import date2num
from iris.util import unify_time_units
import cf_units
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.dates as mdates
import calendar

matplotlib.use('Agg')


def build_cubelist_ensemble(i, array):

    cubes_in_file = iris.load(model_file_paths_ensemble[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        if len(coord_names) == 3:
            array[i] = cube
        if len(coord_names) == 4 and 'depth' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'Level' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'air_pressure' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'height' in coord_names:
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

        """If variable is mrso, distinguish between mrso and mrsos to find model files."""
        if variable == 'mrso':
            variable = 'mrso_'

        """If variable is mrro, distinguish between mrro and mrros to find model files."""
        if variable == 'mrro':
            variable = 'mrro_'

        if variable == 'evaporation':
            variable = 'hfls'

        model_file_paths = find_model_paths_ensemble(list_of_models, model_type, variable, root_directory)

        """If variable is pr_, convert variable back to pr"""
        if variable == 'pr_':
            variable = 'pr'

        """If variable is evspsbl_, convert variable back to evspsbl"""
        if variable == 'evspsbl_':
            variable = 'evspsbl'

        """If variable is mrso, convert variable back to mrso"""
        if variable == 'mrso_':
            variable = 'mrso'

        """If variable is mrro, convert variable back to mrro"""
        if variable == 'mrro_':
            variable = 'mrro'

        if variable == 'hfls':
            variable = 'evaporation'

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

def slicing_ensemble(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        if variable != "treeFrac":
            if variable == 'swc_anom':
                pass
            else:
                cubelist_ensemble[i] = cubelist_ensemble[i].extract(time_range)
            time_points = cubelist_ensemble[i].coord('time').points
            times = cubelist_ensemble[i].coord('time').units.num2date(time_points)

        if variable in ['nrad', 'vpd', 'hurs', 'tas', 'ws', 'sa', 'cancap', 'swc_anom']:
            model_id = cubelist_ensemble[i].long_name
        else:
            model_id = cubelist_ensemble[i].attributes['model_id']
        print model_id
        if variable != "treeFrac":
            print len(times)
            print times[0]
            print times[-1]

    """Slice the regridded cube down to the African domain."""
    cubelist_ensemble[i] = cubelist_ensemble[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

    """Select model ID."""
    if variable in ['nrad', 'vpd', 'hurs', 'tas', 'ws', 'sa', 'cancap', 'swc_anom']:
        model_id = cubelist_ensemble[i].long_name
    else:
        model_id = cubelist_ensemble[i].attributes['model_id']

    if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg' or variable == 'mrros' or variable == 'mrro' or variable == 'prveg':
        cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 86400)

    """if variable is evaporation, divide by 28"""
    if variable == 'evaporation':
        if unit_plot == "mm day-1":
            cubelist_ensemble[i] = iris.analysis.maths.divide(cubelist_ensemble[i], 28)

    """if variable is vpd, multiply by -1"""
    if variable == 'vpd':
        cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], -1)

    if variable == 'evspsblsoi' or variable == "evspsblveg" or variable == 'tran':
        if model_id == "IPSL-CM5B-LR":
            cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 4)
        if unit_plot == 'W m-2':
            cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 28)

    """Reassign model ID."""
    cubelist_ensemble[i].long_name = model_id


    for j in cubelist_ensemble:
        time_points = j.coord('time').points
        times = j.coord('time').units.num2date(time_points)
        #
        # print time_points
        # print "hi"
        # print j

    """Do Hovmoller Calculations"""

    iris.coord_categorisation.add_month_number(cubelist_ensemble[i], 'time', name='month')
    cubelist_ensemble[i] = cubelist_ensemble[i].aggregated_by(['month'], iris.analysis.MEAN)

    array[i] = cubelist_ensemble[i]


def build_cubelist_reanalysis(i, array):

    cubes_in_file = iris.load(reanalysis_file_paths[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        print coord_names
        if len(coord_names) == 3:
            array[i] = cube
        if len(coord_names) == 4 and 'depth' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'Level' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'air_pressure' in coord_names:
            array[i] = cube

def reanalysis_file_paths_func(list_of_reanalysis, variable):
    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles"

    """If variable is pr, distinguish between pr and precipitable water to find model files."""
    if variable == 'pr':
        variable = 'pr_'

    """If variable is evspsbl, distinguish between evspsbl and evspsblsoi to find model files."""
    if variable == 'evspsbl':
        variable = 'evspsbl_'

    """If variable is mrso, distinguish between mrso and mrsos to find model files."""
    if variable == 'mrso':
        variable = 'mrso_'

    """If variable is mrro, distinguish between mrro and mrros to find model files."""
    if variable == 'mrro':
        variable = 'mrro_'

    if variable == 'evaporation':
        variable = 'hfls'

    reanalysis_file_paths = find_reanalysis_file_paths(list_of_reanalysis, variable, root_directory)

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """If variable is evspsbl_, convert variable back to evspsbl"""
    if variable == 'evspsbl_':
        variable = 'evspsbl'

    """If variable is mrso, convert variable back to mrso"""
    if variable == 'mrso_':
        variable = 'mrso'

    """If variable is mrro, convert variable back to mrro"""
    if variable == 'mrro_':
        variable = 'mrro'

    if variable == 'hfls':
        variable = 'evaporation'

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
    print reanalysis_file_paths
    return reanalysis_file_paths


def slicing_reanalysis(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        reanalysis_id = list_of_reanalysis[i]
        print reanalysis_id
        if variable != "treeFrac":
            print cubelist_reanalysis[i]
            if reanalysis_id != "era5":
                if variable == 'swc_anom':
                    pass
                else:
                    cubelist_reanalysis[i] = cubelist_reanalysis[i].extract(time_range)
            if reanalysis_id == "era5":
                cubelist_reanalysis[i] = cubelist_reanalysis[i].extract(iris.Constraint(time=lambda cell: 2010 <= cell.point.year <= 2017))
            time_points = cubelist_reanalysis[i].coord('time').points
            times = cubelist_reanalysis[i].coord('time').units.num2date(time_points)
        if variable != "treeFrac":
            print reanalysis_id
            print len(times)
            print times[0]
            print times[-1]

    """Slice the regridded cube down to the African domain."""
    cubelist_reanalysis[i] = cubelist_reanalysis[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

    """if variable is evaporation, divide by 28"""
    if variable == 'evaporation':
        if unit_plot == 'mm day-1':
            cubelist_reanalysis[i] = iris.analysis.maths.divide(cubelist_reanalysis[i], 28)

    if variable == 'evspsblsoi' or variable == "evspsblveg" or variable == 'tran':
        if unit_plot == 'W m-2':
            cubelist_reanalysis[i] = iris.analysis.maths.multiply(cubelist_reanalysis[i], 28)

    """if variable is vpd, multiply by -1"""
    if variable == 'vpd':
        cubelist_reanalysis[i] = iris.analysis.maths.multiply(cubelist_reanalysis[i], -1)

    """Select the reanalysis ID."""
    reanalysis_id = list_of_reanalysis[i]

    if reanalysis_id == "cfsr":
        cubelist_reanalysis[i].long_name = "CFSR"
        cubelist_reanalysis[i].rename("CFSR")
    if reanalysis_id == "erai":
        cubelist_reanalysis[i].long_name = "ERA-Interim"
        cubelist_reanalysis[i].rename("ERA-Interim")
    if reanalysis_id == 'era5':
        cubelist_reanalysis[i].long_name = "ERA-5"
        cubelist_reanalysis[i].rename("ERA-5")
    if reanalysis_id == "gleam":
        cubelist_reanalysis[i].long_name = "GLEAM"
        cubelist_reanalysis[i].rename("GLEAM")
    if reanalysis_id == "jra":
        cubelist_reanalysis[i].long_name = "JRA-55"
        cubelist_reanalysis[i].rename("JRA-55")
    if reanalysis_id == "merra2":
        cubelist_reanalysis[i].long_name = "MERRA-2"
        cubelist_reanalysis[i].rename("MERRA-2")
    if reanalysis_id == "mswep":
        cubelist_reanalysis[i].long_name = "MSWEP"
        cubelist_reanalysis[i].rename("MSWEP")
    if reanalysis_id == "ncep-ncar":
        cubelist_reanalysis[i].long_name = "NCEP/NCAR"
        cubelist_reanalysis[i].rename("NCEP/NCAR")
    if reanalysis_id == "ncep-doe":
        cubelist_reanalysis[i].long_name = "NCEP DOE-2"
        cubelist_reanalysis[i].rename("NCEP DOE-2")
    if reanalysis_id == 'gewex':
        cubelist_reanalysis[i].long_name = "GEWEX"
        cubelist_reanalysis[i].rename("GEWEX")
    if reanalysis_id == 'modis':
        cubelist_reanalysis[i].long_name = "MODIS"
        cubelist_reanalysis[i].rename("MODIS")
    if reanalysis_id == 'landfluxeval':
        cubelist_reanalysis[i].long_name = "LandFlux-EVAL"
        cubelist_reanalysis[i].rename("LandFlux-EVAL")
    if reanalysis_id == 'chirps':
        cubelist_reanalysis[i].long_name = "CHIRPS"
        cubelist_reanalysis[i].rename("CHIRPS")
    if reanalysis_id == 'trmm':
        cubelist_reanalysis[i].long_name = "TRMM"
        cubelist_reanalysis[i].rename("TRMM")
    if reanalysis_id == 'gpcc':
        cubelist_reanalysis[i].long_name = "GPCC"
        cubelist_reanalysis[i].rename("GPCC")

    print "hi"


    # for j in cubelist_reanalysis:
    #
    #     #
    #     # print time_points
    #     # print "hi"
    #     # print j
    #
    #     """Do Hovmoller Calculations"""
    #
    #     print i

    iris.coord_categorisation.add_month_number(cubelist_reanalysis[i], 'time', name='month')
    cubelist_reanalysis[i] = cubelist_reanalysis[i].aggregated_by(['month'], iris.analysis.MEAN)

    """Format cube to Hovmoller coordinates."""

    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist_reanalysis[i] = cubelist_reanalysis[i] .intersection(longitude=(lower_lon-1, upper_lon+1))

    cubelist_reanalysis[i] = cubelist_reanalysis[i].collapsed('longitude', iris.analysis.MEAN)

    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist_reanalysis[i] = cubelist_reanalysis[i].intersection(latitude=(lower_lat-2, upper_lat+1))


    array[i] = cubelist_reanalysis[i]


def contour_lev_colour_map(lower_value, higher_value, interval, cmap):

    contour_levels = np.arange(lower_value, higher_value+interval, interval)
    cmap = matplotlib.cm.get_cmap(cmap)
    return contour_levels, cmap

def colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, interval):

    print lower_tick
    print upper_tick+interval
    print interval

    if interval >= 1 or interval <= -1:
        colour_bar.set_ticks(np.arange(lower_tick, upper_tick+interval, interval).astype(np.float))
        colour_bar.set_ticklabels(np.arange(lower_tick, upper_tick+interval, interval).astype(np.float))

    if -1 < interval < 1:

        x = np.arange(lower_tick, upper_tick+interval, interval)

        count = 0
        for i in x:
            if -0.0000001 < i < 0.00000001:
                x[count] = 0
            count +=1

        print x

        colour_bar.set_ticks(x)
        colour_bar.set_ticklabels(x)


    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0, labelsize=6)

    return colour_bar

def plot_hovmoller(cube, lower_lat, upper_lat, contour_levels, cmap):

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    contour_plot = iplt.contourf(cube, contour_levels, coords=['month', 'latitude'], cmap=cmap, extend='both')

    #qplt.contourf(ensemble_mean, coords=['month', 'latitude'])

    ax1.set_xlim([1, 12])

    ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
    ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')
    ax1.patch.set_visible(False)

    #contour_plot = iplt.contourf(cube, coords=['month', 'latitude'])

    objects = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec')

    x_pos = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    # print x_pos
    plt.xticks(x_pos, objects, fontsize= 10)

    ax1.set_ylim([lower_lat, upper_lat])

    plt.ylabel('Latitude ($^\circ$)', fontsize = 10)

    colour_bar = plt.colorbar(contour_plot, orientation='horizontal')

    colour_bar = colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, tick_interval)

    if variable == 'tran':
        label = 'Transpiration (mm $\mathregular{day^{-1}}$)'

    if variable == 'nrad':
        label = 'Surface Net Downward Radiation (W $\mathregular{m^{-2}}$)'

    colour_bar.set_label(label, fontsize=9)

    fig.savefig("Hovmoller.png", bbox_inches='tight')



if __name__ == "__main__":

    list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    list_of_models = []
    list_of_reanalysis = ['gewex']

    model_type = "amip"

    list_of_variables = ['nrad']

    # #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008

    lower_value = 80.0
    higher_value = 180.0
    value_interval = 10.0
    lower_tick = 80.0
    upper_tick = 180.0
    tick_interval = 10.0

    cmap = "YlGnBu"

    models = 'no'
    ensemble = 'no'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"

    legend_in_plot = 'yes'
    include_legend = 'yes'


    time_range = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)

    for i in list_of_variables:

        variable = i
        print variable

        if ensemble == 'yes':

            if variable not in ['mrsos', 'lai', 'cancap', 'hurs', 'vpd']:

                """Extract the regridded model file paths for the ensemble mean."""
                model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

                print model_file_paths_ensemble

            """ENSEMBLE ALL MODELS."""

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
                p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            # """Replace fill values with nan."""
            # for i in range(len(cubelist_ensemble)):
            #     array = np.ma.filled(cubelist_ensemble[i].data)
            #     array[array==1e20] = np.nan
            #     cubelist_ensemble[i].data = array
            #     #cubelist_ensemble[i] = cubelist_ensemble[i].data

            """Modify the time units of each cube, to enable ensemble mean to be calculated """

            unify_time_units(cubelist_ensemble)

            print cubelist_ensemble

            for j in cubelist_ensemble:

                j.coord('time').units = cf_units.Unit(j.coord('time').units.origin, calendar='standard')
                time_points = j.coord('time').points
                times = num2date(time_points, units=str(j.coord('time').units))
                times = np.asarray([k.replace(k.year, k.month, 1, 0) for k in times])
                times = date2num(times, units=str(j.coord('time').units))
                j.coord('time').points = times
                j.coord('time').bounds = None
                print j.coord('time')

            """Calculate ensemble mean, and format cube to Hovmoller coordinates."""
            ensemble_mean = sum(cubelist_ensemble) / float(len(cubelist_ensemble))

            with iris.FUTURE.context(cell_datetime_objects=True):
                ensemble_mean = ensemble_mean.intersection(longitude=(lower_lon-1, upper_lon+1))

            ensemble_mean = ensemble_mean.collapsed('longitude', iris.analysis.MEAN)

            with iris.FUTURE.context(cell_datetime_objects=True):
                ensemble_mean = ensemble_mean.intersection(latitude=(lower_lat-2, upper_lat+1))

            single_month = ensemble_mean.extract(iris.Constraint(month=1))

            cube_list = iris.cube.CubeList([ensemble_mean, single_month])

            contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
            cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

            plot_hovmoller(ensemble_mean, lower_lat, upper_lat, contour_levels, cmap)

        if reanalysis == 'yes':

            list_of_reanalysis.sort()

            print list_of_reanalysis

            """Build a list of cubes from the reanalysis file paths."""
            reanalysis_file_paths = reanalysis_file_paths_func(list_of_reanalysis, variable)
            print "hi80"
            print reanalysis_file_paths
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
                p = multiprocessing.Process(target=slicing_reanalysis, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_reanalysis = array.values()
            print "hi2"
            print cubelist_reanalysis[0]

            cube_reanalysis = cubelist_reanalysis[0]

            print cube_reanalysis

            contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
            cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

            plot_hovmoller(cube_reanalysis, lower_lat, upper_lat, contour_levels, cmap)


        # print ensemble_mean
        # print single_month

        #print cube_list[0]

        #print cube_list[1]

        # cube_list.concatenate_cube()
        #
        # print cube_list






        #
        # plot_hovmoller(ensemble_mean)

        #print single_month

        #if duplicate_seasonal_cycle == 'yes':

        # cubelist = []

        #ensemble_mean_2 = ensemble_mean

        #single_month = ensemble_mean.extract(iris.Constraint(month=1))
        #
        # cubelist = np.append(cubelist, ensemble_mean)
        # cubelist = np.append(cubelist, ensemble_mean_2)
        # #cubelist = np.append(cubelist, single_month)

        #cube_list = iris.cube.CubeList([ensemble_mean, ensemble_mean])

        #unify_time_units(cube_list)

        #new_cube = cube_list.concatenate_cube()

        # print new_cube
        #
        # print new_cube.coord('time').points


        # new_cube = cube_list.concatenate()[0]
        #
        # #print cubelist
        #
        # print(cubelist.concatenate())

        #print ensemble_mean_last_month



        # ensemble_mean_last_month.remove_coord('time')
        #
        # print ensemble_mean_last_month
        #
        # ensemble_mean_last_month.add_dim_coord(time_coord, 1)
        #
        # #print ensemble_mean_last_month


        # print cubelist
        #
        # ensemble_mean = cubelist.concatenate()
        #
        # print ensemble_mean
