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
        if reanalysis_id == 'modis' and variable == 'lai':
            pass
        else:
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

    if reanalysis_id == 'modis' and variable == 'lai':
        pass
    else:

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

    if reanalysis_id == 'modis' and variable == 'lai':
        pass
    else:

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

    if interval >= 1 or interval <= -1:

        contour_levels = np.arange(lower_value, higher_value+interval, interval)
        cmap = matplotlib.cm.get_cmap(cmap)
        return contour_levels, cmap

    if -1 < interval < 1:

        start = lower_tick*10
        end = (upper_tick+interval)*10
        new_interval = interval*10

        list = [float(i)/10 for i in np.arange(start, end, new_interval)]

        contour_levels = list
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

        start = lower_tick*10
        end = (upper_tick+interval)*10
        new_interval = interval*10

        list = [float(i)/10 for i in np.arange(start, end, new_interval)]

        print list

        # x = np.arange(lower_tick, upper_tick+interval, interval)

        # count = 0
        # for i in x:
        #     if -0.0000001 < i < 0.00000001:
        #         x[count] = 0
        #     count +=1
        #
        # print x

        colour_bar.set_ticks(list)
        colour_bar.set_ticklabels(list)


    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0, labelsize=6)

    return colour_bar

def plot_hovmoller(cube, lower_lat, upper_lat, contour_levels, cmap, list_of_reanalysis, latitude_points):

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
    ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')
    ax1.patch.set_visible(False)

    if list_of_reanalysis[0] == 'modis' and variable == 'lai':

        print cube.shape
        print latitude_points

        x_pos = np.arange(0, 47, 1)
        objects = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan')

        x,y = np.meshgrid(x_pos, latitude_points)
        #
        contour_plot = plt.contourf(x, y, cube, contour_levels, cmap=cmap, extend='both')

        # x1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
        # squad = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
        #
        #ax1.set_yticklabels(np.linspace(-14, 4, 20))
        # ax1.set_yticklabels(squad, minor=False, rotation=45)

        #
        ax1.set_xticks(np.linspace(1, 46, 13))
        ax1.set_xticklabels(objects)
        ax1.set_xlim([1, 46])
        #
        y_ticks = np.arange(-14.0, 5.0, 2.0)
        ax1.set_yticks(np.linspace(-8.98750019, 42.98749924, 10))
        ax1.set_yticklabels(y_ticks)

    else:

        objects = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan')
        x_pos = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

        x,y = np.meshgrid(x_pos, latitude_points)

        contour_plot = plt.contourf(x, y, cube, contour_levels, cmap=cmap, extend='both')

        print x_pos
        plt.xticks(x_pos, objects, fontsize= 10)

        # ax1.set_yticks(np.linspace(0, 20, 1))
        # print "hi5"
        # print latitude_points
        #ax1.set_yticklabels(latitude_points[::-1])
        ax1.set_ylim([lower_lat, upper_lat])

    plt.ylabel('Latitude ($^\circ$)', fontsize = 10)

    colour_bar = plt.colorbar(contour_plot, orientation='horizontal')

    print lower_tick
    print upper_tick
    print tick_interval

    colour_bar = colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, tick_interval)

    if variable == 'pr':
        label = 'Precipitation (mm $\mathregular{day^{-1}}$)'
    if variable == 'hfls':
        label = 'Surface Upward Latent Heat Flux (W $\mathregular{m^{-2}}$)'
    if variable == 'evapotranspiration':
        label = 'Evapotranspiration (mm $\mathregular{day^{-1}}$)'
    if variable == 'hfss':
        label = 'Surface Upward Sensible Heat Flux (W $\mathregular{m^{-2}}$)'
    if variable == 'evap_fraction':
        label = 'Evaporative Fraction'
    if variable == 'nrad':
        label = 'Surface Net Radiation (W $\mathregular{m^{-2}}$)'
    if variable == 'mrsos':
        label = 'Soil Moisture Content of Upper Layer (mm)'
    if variable == 'mrso':
        label = 'Soil Moisture Content (mm)'
    if variable == 'mrsor':
        label = 'Root-Zone Soil Moisture (m3/m3)'
    if variable == 'tran':
        label = 'Transpiration (mm $\mathregular{day^{-1}}$)'
    if variable == 'evspsbl':
        label = 'Evaporation (mm $\mathregular{day^{-1}}$)'
    if variable == 'evspsblsoi':
        label = 'Bare Soil Evaporation (mm $\mathregular{day^{-1}}$)'
    if variable == 'evspsblveg':
        label = 'Evaporation from Canopy (mm $\mathregular{day^{-1}}$)'
    if variable == 'prveg':
        label = 'Precipitation Intercepted by Canopy (mm $\mathregular{day^{-1}}$)'
    if variable == 'mrros':
        label = 'Surface Runoff Flux (mm $\mathregular{day^{-1}}$)'
    if variable == 'lai':
        label = 'Leaf Area Index'
    if variable == 'mrro':
        label = 'Runoff Flux (mm $\mathregular{day^{-1}}$)'
    if variable == 'rsds':
        label = 'Surface Downward Shortwave Radiation (W $\mathregular{m^{-2}}$)'
    if variable == 'rlds':
        label = 'Surface Downward Longwave Radiation (W $\mathregular{m^{-2}}$)'
    if variable == 'rsus':
        label = 'Surface Upward Shortwave Radiation (W $\mathregular{m^{-2}}$)'
    if variable == 'rlus':
        label = 'Surface Upward Longwave Radiation (W $\mathregular{m^{-2}}$)'
    if variable == 'par':
        label = 'Surface Downward Photosynthetically Active Radiation (W $\mathregular{m^{-2}}$)'

    colour_bar.set_label(label, fontsize=9)

    fig.savefig("Hovmoller.png", bbox_inches='tight')


if __name__ == "__main__":

    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    list_of_models = []
    list_of_reanalysis = ['gewex']
    #list_of_reanalysis = ['none']

    model_type = "amip"

    list_of_variables = ['par']

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

    if list_of_reanalysis[0] == 'modis' and list_of_variables[0] == 'lai':
        pass
    else:
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

                print j.coord('time').units
                j.coord('time').units = cf_units.Unit(j.coord('time').units.origin, calendar='standard')
                time_points = j.coord('time').points
                times = num2date(time_points, units=str(j.coord('time').units))
                times = np.asarray([k.replace(k.year, k.month, 1, 0) for k in times])
                times = date2num(times, units=str(j.coord('time').units))
                j.coord('time').points = times
                j.coord('time').bounds = None
                #print j.coord('time').units

            """Calculate ensemble mean, and format cube to Hovmoller coordinates."""
            ensemble_mean = sum(cubelist_ensemble) / float(len(cubelist_ensemble))

            with iris.FUTURE.context(cell_datetime_objects=True):
                ensemble_mean = ensemble_mean.intersection(longitude=(lower_lon-1, upper_lon+1))

            ensemble_mean = ensemble_mean.collapsed('longitude', iris.analysis.MEAN)

            with iris.FUTURE.context(cell_datetime_objects=True):
                ensemble_mean = ensemble_mean.intersection(latitude=(lower_lat-2, upper_lat+1))

            """Code to add extra month to hovmoller."""
            latitude_points = ensemble_mean.coord('latitude').points[::-1]

            cube_data_all = ensemble_mean.data
            first_month = cube_data_all[0,:]
            len_latitude = len(ensemble_mean.coord('latitude').points)

            ensemble_mean = np.zeros((13,len_latitude))
            ensemble_mean[0:12] = cube_data_all
            ensemble_mean[12] = first_month

            ensemble_mean = ensemble_mean.T[::-1]

            contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
            cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

            plot_hovmoller(ensemble_mean, lower_lat, upper_lat, contour_levels, cmap, list_of_reanalysis, latitude_points)

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
            #print cubelist_reanalysis[0]

            cube_reanalysis = cubelist_reanalysis[0]

            """Code to add extra month to hovmoller."""
            latitude_points = cube_reanalysis.coord('latitude').points[::-1]

            cube_data_all = cube_reanalysis.data
            first_month = cube_data_all[0,:]
            len_latitude = len(cube_reanalysis.coord('latitude').points)

            if list_of_reanalysis[0] == 'modis' and variable == 'lai':
                cube_reanalysis_data = np.zeros((47,len_latitude))
                cube_reanalysis_data[0:46] = cube_data_all
                cube_reanalysis_data[46] = first_month
            else:
                cube_reanalysis_data = np.zeros((13,len_latitude))
                cube_reanalysis_data[0:12] = cube_data_all
                cube_reanalysis_data[12] = first_month

            cube_reanalysis_data = cube_reanalysis_data.T[::-1]

            """Send to plotting function."""
            contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
            cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

            plot_hovmoller(cube_reanalysis_data, lower_lat, upper_lat, contour_levels, cmap, list_of_reanalysis, latitude_points)
