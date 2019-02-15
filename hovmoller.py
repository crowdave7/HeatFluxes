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

        print time_points
        print "hi"
        print j

    """Do Hovmoller Calculations"""

    iris.coord_categorisation.add_month_number(cubelist_ensemble[i], 'time', name='month')
    cubelist_ensemble[i] = cubelist_ensemble[i].aggregated_by(['month'], iris.analysis.MEAN)

    array[i] = cubelist_ensemble[i]

    for j in cubelist_ensemble:
        print "hi2"
        print j.data

    #     time_points = j.coord('time').points
    #     times = j.coord('time').units.num2date(time_points)
    #
    #     print times



if __name__ == "__main__":

    list_of_models = ["bcc-csm1-1-m/"]

    model_type = "amip"

    list_of_variables = ['tran']

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

    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'no'
    plot = 'yes'
    unit_plot = "mm day-1"

    legend_in_plot = 'yes'
    include_legend = 'yes'

    time_range = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)

    for i in list_of_variables:

        variable = i
        print variable

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



            #
            # j.coord('time').units = cf_units.Unit(j.coord('time').units.origin, calendar='standard')
            #
            # print times
            #

            # # print j.coord('time').points
            # #
            # print j.coord('time')

        #
        #
        #     print j.coord('time').units

            #print j.coord('time').units
            #
            # print times
        #
        ensemble_mean = sum(cubelist_ensemble) / float(len(cubelist_ensemble))
        #
        # time_mean = ensemble_mean.collapsed('time', iris.analysis.MEAN)
        #
        # print time_mean
        #
        # qplt.contourf(time_mean)
        #
        # plt.show()


        # with iris.FUTURE.context(cell_datetime_objects=True):
        #     ensemble_mean = ensemble_mean.intersection(longitude=(lower_lon-1, upper_lon+1))
        #
        # print ensemble_mean
        #
        # ensemble_mean = ensemble_mean.collapsed('longitude', iris.analysis.MEAN)
        #
        # with iris.FUTURE.context(cell_datetime_objects=True):
        #     ensemble_mean = ensemble_mean.intersection(latitude=(lower_lat-2, upper_lat+1))
        #
        #
        # print ensemble_mean.data
        #
        # print ensemble_mean
        #
        # qplt.contourf(ensemble_mean, 20)
        #
        # plt.show()
        #


        #print ensemble_mean
