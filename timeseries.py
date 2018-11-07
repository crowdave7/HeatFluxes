#!/usr/bin/env python

"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.analysis.cartography
import iris.util
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans, interp
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
from scipy.stats.stats import pearsonr
import time
matplotlib.use('Agg')

def build_cubelist(i, array):

    cubes_in_file = iris.load(model_file_paths[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        print coord_names
        if len(coord_names) == 3:
            array[i] = cube
        if len(coord_names) == 4 and 'depth' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'air_pressure' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'height' in coord_names:
            array[i] = cube


def model_file_paths_func(list_of_models, model_type, variable):

    if variable == 'nrad':
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1NetRadiationModelFiles"
    if variable == 'vpd' or variable == 'hurs' or variable == 'tas' or variable == 'ws' or variable == 'sa' or variable == 'swc_anom':
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles"
    if variable not in ['nrad', 'vpd', 'hurs', 'tas', 'ws', 'sa', 'swc_anom']:
        root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

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

    """FIND MODEL PATHS."""
    print "finding model paths"
    print root_directory
    model_file_paths = find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory)

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


def find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory):

    if variable in ['nrad', 'vpd', 'hurs', 'tas', 'ws', 'sa', 'cancap', 'swc_anom']:
        print "hi9"
        """Find the paths to the files containing the model data"""
        model_file_paths = []
        print root_directory
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                list_of_models = [i.replace("bcc-csm1-1/", "bcc-csm1-1_") for i in list_of_models]
                list_of_models = [i.replace('/', '') for i in list_of_models]
                for j in list_of_models:
                    if j in path and model_type in path and variable in path:
                        model_file_paths = np.append(model_file_paths, path)

        model_file_paths_nrad_sorted = []

        for i in list_of_models:
            print "hi8"
            print i
            for j in model_file_paths:
                if i in j:
                    print "hi3"
                    print i
                    print j
                    print "hi4"
                    model_file_paths_nrad_sorted = np.append(model_file_paths_nrad_sorted, j)

        print model_file_paths_nrad_sorted
        model_file_paths = model_file_paths_nrad_sorted

    if variable == "evap_fraction":

        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                for j in list_of_models:
                    if j in path and model_type in path and ensemble in path and ('hfss' in path or 'hfls' in path):
                        model_file_paths = np.append(model_file_paths, path)

        model_file_paths_sorted = []

        for i in list_of_models:

            for j in model_file_paths:
                if i in j:
                    print i
                    print j
                    model_file_paths_sorted = np.append(model_file_paths_sorted, j)

    if variable not in ['nrad', 'evap_fraction', 'vpd', 'hurs', 'tas', 'ws', 'sa', 'swc_anom']:

        print list_of_models
        """Find the paths to the directories containing the model data"""
        directory_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in directories:
                path = os.path.join(root, i)
                for j in list_of_models:
                    if j in path and model_type in path and ensemble in path:
                        print path
                        directory_paths = np.append(directory_paths, path)

        """Find the model files and their absolute paths."""
        model_file_paths = []
        for i in directory_paths:
            files = os.listdir(i)
            for j in files:
                if variable in j:
                    model_file_path = os.path.join(i, j)
                    model_file_paths = np.append(model_file_paths, model_file_path)

        #file_names = [os.path.basename(i) for i in model_file_paths]

        #model_file_paths_sorted = sorted(model_file_paths, key=lambda i: os.path.basename(i)[0])

        for i in list_of_models:
            if "CCSM4" in i:
                if variable == "evspsblveg":
                    model_file_paths = np.append(model_file_paths, "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/CCSM4/atlas_evspsblveg_Lmon_CCSM4_amip_r1i1p1_197901-201012.nc")
                if variable == "evspsblsoi":
                    model_file_paths = np.append(model_file_paths, "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/CCSM4/atlas_evspsblsoi_Lmon_CCSM4_amip_r1i1p1_197901-201012.nc")
                if variable == "tran":
                    model_file_paths = np.append(model_file_paths, "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/CCSM4/atlas_tran_Lmon_CCSM4_amip_r1i1p1_197901-201012.nc")

        model_file_paths_sorted = []

        for i in list_of_models:
            for j in model_file_paths:
                if i in j:
                    print i
                    print j
                    model_file_paths_sorted = np.append(model_file_paths_sorted, j)

        count = 0
        for i in model_file_paths_sorted:
            if "GISS-E2-R" in i:
                if "evspsblsoi" in i:
                    model_file_paths_sorted[count] = "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/atlas_evspsblsoi_Lmon_GISS-E2-R_amip_r1i1p1_188001-201012_correct.nc"
                    print model_file_paths_sorted[count]
            count +=1

        count = 0
        for i in model_file_paths_sorted:
            if "BNU-ESM" in i:
                if "evspsblveg" in i:
                    model_file_paths_sorted[count] = '/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/atlas_evspsblveg_Lmon_BNU-ESM_amip_r1i1p1_197901-200812_correct.nc'
                    print model_file_paths_sorted[count]
            count +=1

        count = 0
        for i in model_file_paths_sorted:
            if "CNRM-CM5" in i:
                if "mrsos" in i:
                    model_file_paths_sorted[count] = '/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/atlas_mrsos_Lmon_CNRM-CM5_amip_r1i1p1_197901-200812_correct.nc'
                    print model_file_paths_sorted[count]
            count +=1

        model_file_paths = model_file_paths_sorted

    print "hi2"
    print model_file_paths
    return model_file_paths

def slicing(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        if variable != "treeFrac":
            if variable == 'swc_anom':
                pass
            else:
                cubelist[i] = cubelist[i].extract(time_range)
            time_points = cubelist[i].coord('time').points
            times = cubelist[i].coord('time').units.num2date(time_points)
        if variable in ['nrad', 'vpd', 'hurs', 'tas', 'ws', 'sa', 'cancap', 'swc_anom']:
            model_id = cubelist[i].long_name
        else:
            model_id = cubelist[i].attributes['model_id']
        print model_id
        if variable != "treeFrac":
            print len(times)
            print times[0]
            print times[-1]

    """Remove longitude bounds to allow for intersection."""
    cubelist[i].coord('longitude').bounds = None

    """Slice the regridded cube down to the African domain."""
    cubelist[i] = cubelist[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

    """Select model ID."""
    if variable in ['nrad', 'vpd', 'hurs', 'tas', 'ws', 'sa', 'cancap', 'swc_anom']:
        model_id = cubelist[i].long_name
    else:
        model_id = cubelist[i].attributes['model_id']

    if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg' or variable == 'mrros' or variable == 'mrro' or variable == 'prveg':
        cubelist[i] = iris.analysis.maths.multiply(cubelist[i], 86400)

    """if variable is evaporation, divide by 28"""
    if variable == 'evaporation':
        if unit_plot == 'mm day-1':
            cubelist[i] = iris.analysis.maths.divide(cubelist[i], 28)

    """if variable is vpd, multiply by -1"""
    if variable == 'vpd':
        cubelist[i] = iris.analysis.maths.multiply(cubelist[i], -1)

    print model_id

    """if variable is bare soil evap, canopy evap or transpiration and model is IPSL-CM5B-LR:"""
    if variable == 'evspsblsoi' or variable == "evspsblveg" or variable == 'tran':
        if model_id == "IPSL-CM5B-LR":
            cubelist[i] = iris.analysis.maths.multiply(cubelist[i], 4)
        if unit_plot == 'W m-2':
            cubelist[i] = iris.analysis.maths.multiply(cubelist[i], 28)

    """Reassign model ID."""
    cubelist[i].long_name = model_id
    array[i] = cubelist[i]

def return_cube_transpose_one(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist[i] = cubelist[i].intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        array[i] = cubelist[i]

def return_cube_transpose_two(i, array):
    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist[i] = cubelist[i].intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        array[i] = cubelist[i]

def timeseries_calculations(i, array):

    time_range_1 = iris.Constraint(time=lambda cell: cell.point.month == first_month)

    time_range_2 = iris.Constraint(time=lambda cell: cell.point.month == second_month)

    data_unmasked_1 = first_cube.extract(time_range_1)

    data_unmasked_2 = first_cube.extract(time_range_2)

    time_range_year = iris.Constraint(time=lambda cell: cell.point.year == i)

    data_unmasked_3 = data_unmasked_1.extract(time_range_year)

    data_unmasked_4 = data_unmasked_2.extract(time_range_year)

    #data_unmasked_5 = data_unmasked_3.collapsed('time', iris.analysis.MEAN)

    #data_unmasked_6 = data_unmasked_4.collapsed('time', iris.analysis.MEAN)

    data_unmasked = iris.analysis.maths.subtract(data_unmasked_4, data_unmasked_3)

    """If the first coordinate is longitude,"""
    if coord_names[1] == 'longitude':

        """Set up grid of longitudes and latitudes for basemap."""
        map = Basemap(llcrnrlon=lower_lon, llcrnrlat=lower_lat, urcrnrlon=upper_lon, urcrnrlat=upper_lat, projection='mill')
        longitude = second_cube.coord('longitude').points
        latitude = second_cube.coord('latitude').points
        longitude, latitude = np.meshgrid(longitude, latitude)
        x, y = map(longitude, latitude)

        """Set up grid replacing each gridpoint with a 5x5 grid point."""
        x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*5)
        y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*5)
        x2, y2 = np.meshgrid(x2, y2)

        """Transpose the data to set lat first rather than lon."""
        data_unmasked = np.transpose(data_unmasked.data, (1, 0))

        """Interpolate each grid point of the transposed data into a 5x5 grid. Swap dimensions if wrong way round."""
        try:
            data2 = interp(data_unmasked, x[0], y[:, 0], x2, y2 ,order=1)
        except ValueError:
            data2 = interp(data_unmasked, x[0], np.flipud(y[:, 0]), x2, np.flipud(y2) ,order=1)

        """Mask the oceans on the transposed data."""
        lons2, lats2 = map(x2, y2, inverse=True)
        mdata = maskoceans(lons2, lats2, data2, resolution = 'h', grid = 1.25, inlands=False)

        """Plot figure to check that masking has worked."""
        """
        fig = plt.figure()
        map.drawcoastlines(linewidth=2)
        map.drawcountries(linewidth=2)
        map.drawparallels(np.arange(-50, 60, 10), labels=[1, 0, 0, 0], fontsize=10, linewidth=0.4)
        map.drawmeridians(np.arange(-40, 80, 20), labels=[0, 0, 0, 1], fontsize=10, linewidth=0.4)
        if variable == 'pr':
            contour_levels = np.arange(0, 11, 1)
        if variable == 'hfss':
            contour_levels = np.arange(0, 65, 5)
        if variable == 'hfls':
            contour_levels = np.arange(80, 145, 5)
        contour_plot = map.contourf(x2, y2, mdata, contour_levels, extend='both', cmap = 'YlGnBu')
        colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        model_id = second_cube.long_name
        fig.savefig("mask_"+variable+"_"+str(i)+"_"+model_id+".png")
        print "plot done"

        plt.close()
        """

        """Calculate the mean of the data array (excluding nans in mask) for correlation plot."""
        data = np.nanmean(mdata)

        """Append the data to the array outside the loop to produce the data for the seasonal cycle."""
        array[i] = data

    """If the first coordinate is latitude,"""
    if coord_names[1] == 'latitude':

        """Set up grid of longitudes and latitudes for basemap."""
        map = Basemap(llcrnrlon=lower_lon, llcrnrlat=lower_lat, urcrnrlon=upper_lon, urcrnrlat=upper_lat, projection='mill')
        longitude = data_unmasked.coord('longitude').points
        latitude = data_unmasked.coord('latitude').points
        longitude, latitude = np.meshgrid(longitude, latitude)
        x, y = map(longitude, latitude)

        """Set up grid replacing each gridpoint with a 5x5 grid point."""
        x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*5)
        y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*5)
        x2, y2 = np.meshgrid(x2, y2)

        """Interpolate each grid point of the transposed data into a 5x5 grid. Swap dimensions if wrong way round."""
        try:
            data2 = interp(data_unmasked.data, x[0], y[:, 0], x2, y2 ,order=1)
        except ValueError:
            data2 = interp(data_unmasked.data, x[0], np.flipud(y[:, 0]), x2, np.flipud(y2) ,order=1)

        """Mask the oceans on the transposed data."""
        lons2, lats2 = map(x2, y2, inverse=True)
        mdata = maskoceans(lons2, lats2, data2, resolution = 'h', grid = 1.25, inlands=False)

        """Plot figure to check that masking has worked."""
        """
        fig = plt.figure()
        map.drawcoastlines(linewidth=2)
        map.drawcountries(linewidth=2)
        map.drawparallels(np.arange(-50, 60, 10), labels=[1, 0, 0, 0], fontsize=10, linewidth=0.4)
        map.drawmeridians(np.arange(-40, 80, 20), labels=[0, 0, 0, 1], fontsize=10, linewidth=0.4)
        if variable == 'pr':
            contour_levels = np.arange(0, 11, 1)
        if variable == 'hfss':
            contour_levels = np.arange(0, 65, 5)
        if variable == 'hfls':
            contour_levels = np.arange(80, 145, 5)
        contour_plot = map.contourf(x2, y2, mdata, contour_levels, extend='both', cmap = 'YlGnBu')
        colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        model_id = second_cube.long_name
        fig.savefig("mask_"+variable+"_"+str(i)+"_"+model_id+".png")
        print "plot done"

        plt.close()#
        """

        """Calculate the mean of the data array (excluding nans in mask) for correlation plot."""
        data = np.nanmean(mdata)

        """Append the data to the array outside the loop to produce the data for the seasonal cycle."""
        array[i] = data

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

    array[i] = cubelist_ensemble[i]

def return_cube_transpose_one_ensemble(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist_ensemble[i] = cubelist_ensemble[i].intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        array[i] = cubelist_ensemble[i]

def return_cube_transpose_two_ensemble(i, array):
    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist_ensemble[i] = cubelist_ensemble[i].intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        array[i] = cubelist_ensemble[i]

def build_cubelist_reanalysis(i, array):

    cubes_in_file = iris.load(reanalysis_file_paths[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        print cube
        coord_names = [coord.name() for coord in cube.coords()]
        print coord_names
        if len(coord_names) == 3:
            array[i] = cube
        if len(coord_names) == 4 and 'depth' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'Level' in coord_names:
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



    array[i] = cubelist_reanalysis[i]


def return_cube_transpose_one_reanalysis(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist_reanalysis[i] = cubelist_reanalysis[i].intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        array[i] = cubelist_reanalysis[i]

def return_cube_transpose_two_reanalysis(i, array):
    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist_reanalysis[i] = cubelist_reanalysis[i].intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        array[i] = cubelist_reanalysis[i]


def line_colours(number_of_models, number_of_reanalysis, cmap):

    cmap = plt.get_cmap(cmap)
    model_line_colours = [cmap(i) for i in np.linspace(0, 1, number_of_models)]
    reanalysis_line_colours = [cmap(i) for i in np.linspace(0, 1, number_of_reanalysis)]
    return model_line_colours, reanalysis_line_colours

if __name__ == "__main__":

    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/"]
    list_of_reanalysis = []

    model_type = "amip"

    list_of_variables = ['nrad', 'tran']

    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ['BNU-ESM/', 'MIROC5/', "MRI-AGCM3-2S/"]

    #CONGO BASIN
    # lower_lat = -14
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    lower_lat = -14
    upper_lat = -5
    lower_lon = 18
    upper_lon = 29

    lower_year = 1979
    upper_year = 2008

    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'no'
    unit_plot = "mm day-1"

    first_month = 7
    second_month = 11

    time_range = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)

    timeseries_models_variable1_array = []
    timeseries_models_variable2_array = []
    timeseries_ensemble_variable1_array = []
    timeseries_ensemble_variable2_array = []
    timeseries_dv_variable1_array = []
    timeseries_dv_variable2_array = []
    timeseries_pv_variable1_array = []
    timeseries_pv_variable2_array = []
    timeseries_reanalysis_variable1_array = []
    timeseries_reanalysis_variable2_array = []

    count = 0
    for i in list_of_variables:

        timeseries_models_array = []
        timeseries_ensemble_array = []
        timeseries_dv_array = []
        timeseries_pv_array = []
        timeseries_reanalysis_array = []
        model_strings_for_plot = []
        ensemble_string_for_plot = []
        reanalysis_strings_for_plot = []

        variable = i
        print variable

        # ------------------------------------------------------------------------------------------------------------------------------------------

        """MODELS"""

        if len(list_of_models) == 0:
            pass
        else:

            if models == 'no':
                pass
            else:

                """Extract the model file paths."""
                start_time = time.time()
                model_file_paths = model_file_paths_func(list_of_models, model_type, variable)
                print time.time() - start_time, "seconds"

                print len(model_file_paths)

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

                unit = cubelist[0].units


                """Slice each cube by the time range."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(cubelist)):
                    p = multiprocessing.Process(target=slicing, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist = array.values()

                """Define first set of cubes to transpose."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(cubelist)):
                    p = multiprocessing.Process(target=return_cube_transpose_one, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                first_set_transposed_cubes = array.values()

                """Define second set of cubes to transpose."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(cubelist)):
                    p = multiprocessing.Process(target=return_cube_transpose_two, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                second_set_transposed_cubes = array.values()

                """For each cube,"""
                for i in np.arange(0, len(cubelist)):

                    model_id = cubelist[i].long_name

                    with iris.FUTURE.context(cell_datetime_objects=True):

                        """Select the cubes to transpose and the coordinates of the first cube."""
                        first_cube = first_set_transposed_cubes[i]
                        second_cube = second_set_transposed_cubes[i]
                        coord_names = [coord.name() for coord in first_cube.coords()]
                        print "hi8"

                        """Multiprocess calculations for each month."""
                        manager = multiprocessing.Manager()
                        array = manager.dict()
                        jobs = []
                        for i in np.arange(lower_year, upper_year+1, 1):
                            p = multiprocessing.Process(target=timeseries_calculations, args=(i, array))
                            jobs.append(p)
                            p.start()
                        for process in jobs:
                            process.join()
                        model_timeseries = array.values()
                        print model_id
                        model_strings_for_plot = np.append(model_strings_for_plot, model_id)
                        # first_value = model_timeseries[0]
                        # model_timeseries.append(first_value)
                        timeseries_models_array.append(model_timeseries)

                if len(list_of_models) > 0:

                    if models == 'no':
                        pass

                    else:

                        for i in np.arange(0, len(list_of_models)):

                            model_number = i

                            model_id = model_strings_for_plot[i]
                            print model_id
                            model_timeseries = timeseries_models_array[i]

                            model_timeseries = [float('%.3f' % i) for i in model_timeseries]

                            print model_timeseries

                            mean = np.mean(model_timeseries)
                            print mean

                            if count == 0:
                                timeseries_models_variable1_array = timeseries_models_array

                            if count == 1:
                                print "hi16"
                                timeseries_models_variable2_array = timeseries_models_array

        """Ensemble MEAN"""

        if ensemble == 'yes':

            model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

            model_file_paths_ensemble_dv = model_file_paths_ensemble_func(list_dynamic_veg_models, model_type, variable)

            model_file_paths_ensemble_pv = model_file_paths_ensemble_func(list_prescribed_veg_models, model_type, variable)

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

            """Replace fill values with nan."""
            for i in range(len(cubelist_ensemble)):
                array = np.ma.filled(cubelist_ensemble[i].data)
                array[array==1e20] = np.nan
                cubelist_ensemble[i].data = array
                #cubelist_ensemble[i] = cubelist_ensemble[i].data

            """Define first set of cubes to transpose."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=return_cube_transpose_one_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            first_set_transposed_cubes_ensemble = array.values()

            """Define second set of cubes to transpose."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=return_cube_transpose_two_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            second_set_transposed_cubes_ensemble = array.values()

            """For each cube,"""
            for i in np.arange(0, len(cubelist_ensemble)):

                model_id = cubelist_ensemble[i].long_name
                print model_id

                with iris.FUTURE.context(cell_datetime_objects=True):

                    """Select the cubes to transpose and the coordinates of the first cube."""
                    first_cube = first_set_transposed_cubes_ensemble[i]
                    second_cube = second_set_transposed_cubes_ensemble[i]
                    coord_names = [coord.name() for coord in first_cube.coords()]

                    """Multiprocess calculations for each month."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(lower_year, upper_year+1, 1):
                        p = multiprocessing.Process(target=timeseries_calculations, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    ensemble_timeseries = array.values()
                    timeseries_ensemble_array.append(ensemble_timeseries)

                    if count == 0:
                        timeseries_ensemble_variable1_array = timeseries_ensemble_array
                        timeseries_ensemble_variable1_array = np.mean(timeseries_ensemble_variable1_array, axis = 0)

                    if count == 1:
                        timeseries_ensemble_variable2_array = timeseries_ensemble_array
                        timeseries_ensemble_variable2_array = np.mean(timeseries_ensemble_variable2_array, axis = 0)


                ensemble_string_for_plot = np.append(ensemble_string_for_plot, "Ensemble")

        if include_dynamic_prescribed == 'yes':

            model_file_paths_ensemble = model_file_paths_ensemble_dv

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

            print cubelist_ensemble

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

            """Replace fill values with nan."""
            for i in range(len(cubelist_ensemble)):
                array = np.ma.filled(cubelist_ensemble[i].data)
                array[array==1e20] = np.nan
                cubelist_ensemble[i].data = array
                #cubelist_ensemble[i] = cubelist_ensemble[i].data

            """Define first set of cubes to transpose."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=return_cube_transpose_one_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            first_set_transposed_cubes_ensemble = array.values()

            """Define second set of cubes to transpose."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=return_cube_transpose_two_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            second_set_transposed_cubes_ensemble = array.values()

            """For each cube,"""
            for i in np.arange(0, len(cubelist_ensemble)):

                model_id = cubelist_ensemble[i].long_name

                with iris.FUTURE.context(cell_datetime_objects=True):

                    """Select the cubes to transpose and the coordinates of the first cube."""
                    first_cube = first_set_transposed_cubes_ensemble[i]
                    second_cube = second_set_transposed_cubes_ensemble[i]
                    coord_names = [coord.name() for coord in first_cube.coords()]

                    """Multiprocess calculations for each month."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(lower_year, upper_year+1, 1):
                        p = multiprocessing.Process(target=timeseries_calculations, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    dv_timeseries = array.values()
                    timeseries_dv_array.append(dv_timeseries)

                if count == 0:
                    timeseries_dv_variable1_array = timeseries_dv_array
                    timeseries_dv_variable1_array = np.mean(timeseries_dv_variable1_array, axis = 0)

                if count == 1:
                    print "hi16"
                    timeseries_dv_variable2_array = timeseries_dv_array
                    timeseries_dv_variable2_array = np.mean(timeseries_dv_variable2_array, axis = 0)


            model_file_paths_ensemble = model_file_paths_ensemble_pv

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

            print cubelist_ensemble

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

            """Replace fill values with nan."""
            for i in range(len(cubelist_ensemble)):
                array = np.ma.filled(cubelist_ensemble[i].data)
                array[array==1e20] = np.nan
                cubelist_ensemble[i].data = array
                #cubelist_ensemble[i] = cubelist_ensemble[i].data

            """Define first set of cubes to transpose."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=return_cube_transpose_one_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            first_set_transposed_cubes_ensemble = array.values()

            """Define second set of cubes to transpose."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=return_cube_transpose_two_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            second_set_transposed_cubes_ensemble = array.values()

            """For each cube,"""
            for i in np.arange(0, len(cubelist_ensemble)):

                model_id = cubelist_ensemble[i].long_name

                with iris.FUTURE.context(cell_datetime_objects=True):

                    """Select the cubes to transpose and the coordinates of the first cube."""
                    first_cube = first_set_transposed_cubes_ensemble[i]
                    second_cube = second_set_transposed_cubes_ensemble[i]
                    coord_names = [coord.name() for coord in first_cube.coords()]

                    """Multiprocess calculations for each month."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(lower_year, upper_year+1, 1):
                        p = multiprocessing.Process(target=timeseries_calculations, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    pv_timeseries = array.values()
                    timeseries_pv_array.append(pv_timeseries)

                if count == 0:
                    timeseries_pv_variable1_array = timeseries_pv_array
                    timeseries_pv_variable1_array = np.mean(timeseries_pv_variable1_array, axis = 0)

                if count == 1:
                    print "hi16"
                    timeseries_pv_variable2_array = timeseries_pv_array
                    timeseries_pv_variable2_array = np.mean(timeseries_pv_variable2_array, axis = 0)

        """REANALYSES."""

        if len(list_of_reanalysis) == 0:
            pass

        else:

            list_of_reanalysis.sort()

            print list_of_reanalysis

            """Build a list of cubes from the reanalysis file paths."""
            reanalysis_file_paths = reanalysis_file_paths_func(list_of_reanalysis, variable)
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

            """Define first set of cubes to transpose."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_reanalysis)):
                p = multiprocessing.Process(target=return_cube_transpose_one_reanalysis, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            first_set_transposed_cubes_reanalysis = array.values()

            """Define second set of cubes to transpose."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_reanalysis)):
                p = multiprocessing.Process(target=return_cube_transpose_two_reanalysis, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            second_set_transposed_cubes_reanalysis = array.values()

            """For each cube,"""
            for i in np.arange(0, len(cubelist_reanalysis)):

                reanalysis_id = cubelist_reanalysis[i].long_name

                with iris.FUTURE.context(cell_datetime_objects=True):

                    """Select the cubes to transpose and the coordinates of the first cube."""
                    first_cube = first_set_transposed_cubes_reanalysis[i]
                    second_cube = second_set_transposed_cubes_reanalysis[i]
                    coord_names = [coord.name() for coord in first_cube.coords()]

                    """Multiprocess calculations for each month."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    if reanalysis_id == "gleam" and variable == 'tran':
                        lower_year = 1984
                    for i in np.arange(lower_year, upper_year+1, 1):
                        p = multiprocessing.Process(target=timeseries_calculations, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    reanalysis_timeseries = array.values()
                    print reanalysis_timeseries
                    print reanalysis_id
                    reanalysis_strings_for_plot = np.append(reanalysis_strings_for_plot, reanalysis_id)
                    # first_value = model_timeseries[0]
                    # model_timeseries.append(first_value)
                    timeseries_reanalysis_array.append(reanalysis_timeseries)

                if len(list_of_reanalysis) > 0:

                    if reanalysis == 'no':
                        pass

                    else:

                        for i in np.arange(0, len(list_of_reanalysis)):

                            model_number = i

                            reanalysis_id = reanalysis_strings_for_plot[i]
                            print reanalysis_id
                            reanalysis_timeseries = timeseries_reanalysis_array[i]

                            reanalysis_timeseries = [float('%.3f' % i) for i in reanalysis_timeseries]

                            print reanalysis_timeseries

                            mean = np.mean(reanalysis_timeseries)
                            print mean

                            if count == 0:
                                timeseries_reanalysis_variable1_array = timeseries_reanalysis_array

                            if count == 1:
                                print "hi16"
                                timeseries_reanalysis_variable2_array = timeseries_reanalysis_array

        count +=1

    list_all_datasets = []

    if len(list_of_models) > 0:
        list_all_datasets = np.append(list_all_datasets, list_of_models)

    if ensemble == "yes":
        list_all_datasets = np.append(list_all_datasets, "Ensemble")

    if include_dynamic_prescribed == 'yes':
        list_all_datasets = np.append(list_all_datasets, "Best")
        list_all_datasets = np.append(list_all_datasets, "Worst")

    if len(list_of_reanalysis) > 0:
        list_all_datasets = np.append(list_all_datasets, list_of_reanalysis)

    list_correlation_coeff = []

    if len(list_of_models) > 0:

        if models == 'no':
            pass

        else:

            for i in np.arange(0, len(list_of_models)):

                model_number = i
                model_id = model_strings_for_plot[i]
                print model_id
                timeseries_1 = timeseries_models_variable1_array[i]
                mean_1 = np.mean(timeseries_1)
                timeseries_2 = timeseries_models_variable2_array[i]
                mean_2 = np.mean(timeseries_2)
                # print timeseries_1
                print mean_1
                # print timeseries_2
                print mean_2

                pearson = pearsonr(timeseries_1, timeseries_2)
                pearsoncoeff = pearson[0]
                print pearsoncoeff
                list_correlation_coeff = np.append(list_correlation_coeff, pearsoncoeff)

    if ensemble == 'yes':

        print "ensemble"
        timeseries_1 = timeseries_ensemble_variable1_array
        mean_1 = np.mean(timeseries_1)
        print mean_1
        timeseries_2 = timeseries_ensemble_variable2_array
        mean_2 = np.mean(timeseries_2)
        print mean_2

        pearson = pearsonr(timeseries_1, timeseries_2)
        pearsoncoeff = pearson[0]
        print pearsoncoeff
        list_correlation_coeff = np.append(list_correlation_coeff, pearsoncoeff)

    if include_dynamic_prescribed == 'yes':

        print "Best"
        timeseries_1 = timeseries_dv_variable1_array
        mean_1 = np.mean(timeseries_1)
        print mean_1
        timeseries_2 = timeseries_dv_variable2_array
        mean_2 = np.mean(timeseries_2)
        print mean_2

        pearson = pearsonr(timeseries_1, timeseries_2)
        pearsoncoeff = pearson[0]
        print pearsoncoeff
        list_correlation_coeff = np.append(list_correlation_coeff, pearsoncoeff)

        print "Worst"
        timeseries_1 = timeseries_pv_variable1_array
        mean_1 = np.mean(timeseries_1)
        print mean_1
        timeseries_2 = timeseries_pv_variable2_array
        mean_2 = np.mean(timeseries_2)
        print mean_2

        pearson = pearsonr(timeseries_1, timeseries_2)
        pearsoncoeff = pearson[0]
        print pearsoncoeff
        list_correlation_coeff = np.append(list_correlation_coeff, pearsoncoeff)

    if len(list_of_reanalysis) > 0:

        if reanalysis == 'no':
            pass

        else:

            for i in np.arange(0, len(list_of_reanalysis)):

                model_number = i
                reanalyis_id = reanalysis_strings_for_plot[i]
                print reanalysis_id
                timeseries_1 = timeseries_reanalysis_variable1_array[i]
                mean_1 = np.mean(timeseries_1)
                timeseries_2 = timeseries_reanalysis_variable2_array[i]
                mean_2 = np.mean(timeseries_2)
                print timeseries_1
                print mean_1
                print timeseries_2
                print mean_2

                if len(timeseries_1) == 24 and len(timeseries_2) == 29:
                    timeseries_2 = timeseries_2[5:]
                    print len(timeseries_2)
                    print timeseries_2

                pearson = pearsonr(timeseries_1, timeseries_2)
                pearsoncoeff = pearson[0]
                print pearsoncoeff
                list_correlation_coeff = np.append(list_correlation_coeff, pearsoncoeff)


    print list_all_datasets
    print list_correlation_coeff
