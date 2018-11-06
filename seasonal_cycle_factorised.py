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
    print cubelist[i]
    array[i] = cubelist[i]

def return_cube_transpose_one(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist[i] = cubelist[i].intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        array[i] = cubelist[i]

def return_cube_transpose_two(i, array):
    with iris.FUTURE.context(cell_datetime_objects=True):
        cubelist[i] = cubelist[i].intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        array[i] = cubelist[i]

def seasonal_cycle_calculations(i, array):
    """Constrain the data for each month only."""
    time_range = iris.Constraint(time=lambda cell: cell.point.month == i)
    data_unmasked = first_cube.extract(time_range)
    print data_unmasked
    print variable
    print reanalysis_id
    if variable == 'swc_anom':
        pass
    # if variable == 'lai' and reanalysis_id == 'MODIS':
    #     print "hi"
    #     pass
    else:
        data_unmasked = data_unmasked.collapsed('time', iris.analysis.MEAN)

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


def plot_bar(barchart_array, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, bar_season, bar_width, variables_to_subtract):

    if number_of_variables == 3:

        fig = plt.figure()

        variable_data_1 = barchart_array[0,:]
        variable_data_2 = barchart_array[1,:]
        variable_data_3 = barchart_array[2,:]

        if ensemble == 'yes':
            number_of_models = number_of_models+1

            if include_dynamic_prescribed == 'yes':
                number_of_models = number_of_models+2

            error_bar_array_1 = np.zeros(number_of_models)
            standard_dev_var_1 = np.std(variable_data_1)
            error_bar_array_1[-1] = standard_dev_var_1

            error_bar_array_2 = np.zeros(number_of_models)
            standard_dev_var_2 = np.std(variable_data_2)
            error_bar_array_2[-1] = standard_dev_var_2

            error_bar_array_3 = np.zeros(number_of_models)
            standard_dev_var_3 = np.std(variable_data_3)
            error_bar_array_3[-1] = standard_dev_var_3

        if ensemble != 'yes':
            error_bar_array_1 = np.zeros(number_of_models)
            error_bar_array_2 = np.zeros(number_of_models)
            error_bar_array_3 = np.zeros(number_of_models)

        print number_of_models
        print number_of_reanalysis

        number_of_bars = number_of_models + number_of_reanalysis

        number_of_bars = np.arange(number_of_bars)

        print number_of_bars
        print variable_data_1
        print bar_width

        print len(variable_data_1)

        p1 = plt.bar(number_of_bars, variable_data_1, bar_width, color=bar_colours[0])

        for bar in p1:
            bar.set_edgecolor("black")
            bar.set_linewidth(0.5)
            bar.set_linestyle("-")

        #plt.errorbar(number_of_bars-0.07, variable_data_1, yerr=error_bar_array_1, fmt='none', ecolor='black')

        p2 = plt.bar(number_of_bars, variable_data_2, bar_width, bottom=variable_data_1, color=bar_colours[1])

        for bar in p2:
            bar.set_edgecolor("black")
            bar.set_linewidth(0.5)
            bar.set_linestyle("-")

        #plt.errorbar(number_of_bars, variable_data_1 + variable_data_2, yerr=error_bar_array_2, fmt='none', ecolor='black')

        p3 = plt.bar(number_of_bars, variable_data_3, bar_width, hatch="///")

        matplotlib.rcParams['hatch.linewidth'] = 0.4

        for bar in p3:
            bar.set_edgecolor("black")
            bar.set_facecolor("None")
            bar.set_linewidth(0.5)
            bar.set_linestyle("-")

        if include_dynamic_prescribed == 'yes':
            model_strings_for_plot = np.append(model_strings_for_plot, "Strong GCMs")
            model_strings_for_plot = np.append(model_strings_for_plot, "Weak GCMs")

        if ensemble == 'yes':
            model_strings_for_plot = np.append(model_strings_for_plot, "Ensemble")
        if number_of_reanalysis > 0:
            model_strings_for_plot = np.append(model_strings_for_plot, reanalysis_strings_for_plot)

        plt.ylabel(bar_season+" "+bar_y_axis_title, fontsize=9)
        plt.xticks(number_of_bars, model_strings_for_plot, rotation=45, ha='right', fontsize=8)

        plt.ylim(lower_y_lim, upper_y_lim)

        for i in list_of_variables:
            if i == 'tran':
                list_of_variables = [i.replace('tran', 'Transpiration') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'evspsblsoi':
                list_of_variables = [i.replace('evspsblsoi', 'Bare Soil Evaporation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'evspsblveg':
                list_of_variables = [i.replace('evspsblveg', 'Canopy Evaporation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'evaporation':
                list_of_variables = [i.replace('evaporation', 'Evaporation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'nrad':
                list_of_variables = [i.replace('nrad', 'Surface Net Downward Radiation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'hfls':
                list_of_variables = [i.replace('hfls', 'Surface Upward Latent Heat Flux') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'hfss':
                list_of_variables = [i.replace('hfss', 'Surface Upward Sensible Heat Flux') for i in list_of_variables]

        plt.legend((p3[0], p2[0], p1[0]), list_of_variables[::-1], fontsize=9)

        plt.gcf().subplots_adjust(bottom=0.3)

        if bar_season == 'DJF':
            bar_season = 'S1'
        if bar_season == 'MAM':
            bar_season = 'S2'
        if bar_season == 'JJA':
            bar_season = 'S3'
        if bar_season == 'SON':
            bar_season = 'S4'

        print "Saving figure"
        fig.savefig("Bar_Graph_"+bar_season, bbox_inches='tight')


    if number_of_variables == 4:

        print number_of_variables

        fig = plt.figure()

        variable_data_1 = barchart_array[0,:]
        variable_data_2 = barchart_array[1,:]
        variable_data_3 = barchart_array[2,:]
        variable_data_4 = barchart_array[3,:]

        if ensemble == 'yes':
            number_of_models = number_of_models+1

            if include_dynamic_prescribed == 'yes':
                number_of_models = number_of_models+2

            error_bar_array_1 = np.zeros(number_of_models)
            standard_dev_var_1 = np.std(variable_data_1)
            error_bar_array_1[-1] = standard_dev_var_1

            error_bar_array_2 = np.zeros(number_of_models)
            standard_dev_var_2 = np.std(variable_data_2)
            error_bar_array_2[-1] = standard_dev_var_2

            error_bar_array_3 = np.zeros(number_of_models)
            standard_dev_var_3 = np.std(variable_data_3)
            error_bar_array_3[-1] = standard_dev_var_3

        if ensemble != 'yes':
            error_bar_array_1 = np.zeros(number_of_models)
            error_bar_array_2 = np.zeros(number_of_models)
            error_bar_array_3 = np.zeros(number_of_models)

        number_of_bars = number_of_models + number_of_reanalysis

        number_of_bars = np.arange(number_of_bars)

        print variable_data_1
        print variable_data_2
        print variable_data_3
        print variable_data_4

        p1 = plt.bar(number_of_bars, variable_data_1, bar_width, color=bar_colours[0])

        for bar in p1:
            bar.set_edgecolor("black")
            bar.set_linewidth(0.5)
            bar.set_linestyle("-")

        p2 = plt.bar(number_of_bars, variable_data_2, bar_width, bottom=variable_data_1, color=bar_colours[1])

        for bar in p2:
            bar.set_edgecolor("black")
            bar.set_linewidth(0.5)
            bar.set_linestyle("-")

        #plt.errorbar(number_of_bars-0.07, variable_data_1, yerr=error_bar_array_1, fmt='none', ecolor='black')

        p3 = plt.bar(number_of_bars, variable_data_3, bar_width, bottom=np.array(variable_data_1)+np.array(variable_data_2), color=bar_colours[2])

        for bar in p3:
            bar.set_edgecolor("black")
            bar.set_linewidth(0.5)
            bar.set_linestyle("-")

        #plt.errorbar(number_of_bars, variable_data_1 + variable_data_2, yerr=error_bar_array_2, fmt='none', ecolor='black')

        p4 = plt.bar(number_of_bars, variable_data_4, bar_width, hatch="///", zorder = 2)

        matplotlib.rcParams['hatch.linewidth'] = 0.4

        for bar in p4:
            bar.set_edgecolor("black")
            bar.set_facecolor("None")
            bar.set_linewidth(0.5)
            bar.set_linestyle("-")

        #plt.errorbar(number_of_bars+0.13, variable_data_1 + variable_data_2 + variable_data_3, yerr=error_bar_array_3, fmt='none', ecolor='black')

        if len(variables_to_subtract) > 0:

            print barchart_array

            for i in variables_to_subtract:

                if i == ['pr', 'evspsblsoi', 'evspsblveg', 'mrro']:
                    string = 'Available Water for Plant Uptake'
                    color = 'blue'

                if i == ['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']:
                    string = 'Residual Evaporation'
                    color = 'red'

                indices = []
                for j in i:
                    print list_of_variables.index(j)
                    indices = np.append(indices, list_of_variables.index(j))

                print indices

                barchart_array_to_subtract = np.zeros((len(i), len(number_of_bars)))

                print barchart_array_to_subtract

                count = 0

                for k in indices:
                    barchart_array_to_subtract[count] = barchart_array[int(k)]
                    count +=1

                print barchart_array_to_subtract
                print barchart_array_to_subtract.shape
            #
                print len(number_of_bars)
                bars_subtracted = np.zeros((len(number_of_bars)))

                print bars_subtracted.shape

                for i in number_of_bars:

                    bar_to_subtract = barchart_array_to_subtract[:, int(i)]

                    print bar_to_subtract

                    b1 = np.subtract(bar_to_subtract[0], bar_to_subtract[1])
                    b2 = np.subtract(b1, bar_to_subtract[2])
                    bar_subtracted = np.subtract(b2, bar_to_subtract[3])

                    bars_subtracted[i] = bar_subtracted

                print bars_subtracted

                print variable_data_1
                print variable_data_2
                print variable_data_3
                #print color

                p5 = plt.bar(number_of_bars, bars_subtracted, bar_width, bottom=np.array(variable_data_1)+np.array(variable_data_2)+np.array(variable_data_3), color=color)

                for bar in p5:
                    bar.set_edgecolor("black")
                    bar.set_linewidth(0.5)
                    bar.set_linestyle("-")

    if number_of_variables == 5:

        print number_of_variables

        if divide_four_variables_by_fifth == 'yes':

            fig = plt.figure()

            variable_data_1 = barchart_array[0,:]
            variable_data_2 = barchart_array[1,:]
            variable_data_3 = barchart_array[2,:]
            variable_data_4 = barchart_array[3,:]
            variable_data_5 = barchart_array[4,:]

            if ensemble == 'yes':
                number_of_models = number_of_models+1

                if include_dynamic_prescribed == 'yes':
                    number_of_models = number_of_models+2

                error_bar_array_1 = np.zeros(number_of_models)
                standard_dev_var_1 = np.std(variable_data_1)
                error_bar_array_1[-1] = standard_dev_var_1

                error_bar_array_2 = np.zeros(number_of_models)
                standard_dev_var_2 = np.std(variable_data_2)
                error_bar_array_2[-1] = standard_dev_var_2

                error_bar_array_3 = np.zeros(number_of_models)
                standard_dev_var_3 = np.std(variable_data_3)
                error_bar_array_3[-1] = standard_dev_var_3

            if ensemble != 'yes':
                error_bar_array_1 = np.zeros(number_of_models)
                error_bar_array_2 = np.zeros(number_of_models)
                error_bar_array_3 = np.zeros(number_of_models)

            number_of_bars = number_of_models + number_of_reanalysis

            number_of_bars = np.arange(number_of_bars)

            print variable_data_1
            print variable_data_2
            print variable_data_3
            print variable_data_4

            variable_data_1 = variable_data_1 / variable_data_5

            variable_data_2 = variable_data_2 / variable_data_5

            variable_data_3 = variable_data_3 / variable_data_5

            variable_data_4 = variable_data_4 / variable_data_5

            p1 = plt.bar(number_of_bars, variable_data_1, bar_width, color=bar_colours[0])

            for bar in p1:
                bar.set_edgecolor("black")
                bar.set_linewidth(0.5)
                bar.set_linestyle("-")

            p2 = plt.bar(number_of_bars, variable_data_2, bar_width, bottom=variable_data_1, color=bar_colours[1])

            for bar in p2:
                bar.set_edgecolor("black")
                bar.set_linewidth(0.5)
                bar.set_linestyle("-")

            #plt.errorbar(number_of_bars-0.07, variable_data_1, yerr=error_bar_array_1, fmt='none', ecolor='black')

            p3 = plt.bar(number_of_bars, variable_data_3, bar_width, bottom=np.array(variable_data_1)+np.array(variable_data_2), color=bar_colours[2])

            for bar in p3:
                bar.set_edgecolor("black")
                bar.set_linewidth(0.5)
                bar.set_linestyle("-")

            #plt.errorbar(number_of_bars, variable_data_1 + variable_data_2, yerr=error_bar_array_2, fmt='none', ecolor='black')

            p4 = plt.bar(number_of_bars, variable_data_4, bar_width, hatch="///", zorder = 2)

            matplotlib.rcParams['hatch.linewidth'] = 0.4

            for bar in p4:
                bar.set_edgecolor("black")
                bar.set_facecolor("None")
                bar.set_linewidth(0.5)
                bar.set_linestyle("-")

            #plt.errorbar(number_of_bars+0.13, variable_data_1 + variable_data_2 + variable_data_3, yerr=error_bar_array_3, fmt='none', ecolor='black')

            if len(variables_to_subtract) > 0:

                print barchart_array

                for i in variables_to_subtract:

                    if i == ['pr', 'evspsblsoi', 'evspsblveg', 'mrro']:
                        string = 'Available Water for Plant Uptake'
                        color = 'blue'

                    if i == ['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']:
                        string = 'Residual Evaporation'
                        color = 'red'

                    indices = []
                    for j in i:
                        print list_of_variables.index(j)
                        indices = np.append(indices, list_of_variables.index(j))

                    print indices

                    barchart_array_to_subtract = np.zeros((len(i), len(number_of_bars)))

                    print barchart_array_to_subtract

                    count = 0

                    for k in indices:
                        barchart_array_to_subtract[count] = barchart_array[int(k)]
                        count +=1

                    print barchart_array_to_subtract
                    print barchart_array_to_subtract.shape
                #
                    print len(number_of_bars)
                    bars_subtracted = np.zeros((len(number_of_bars)))

                    print bars_subtracted.shape

                    for i in number_of_bars:

                        bar_to_subtract = barchart_array_to_subtract[:, int(i)]

                        print bar_to_subtract

                        b1 = np.subtract(bar_to_subtract[0], bar_to_subtract[1])
                        b2 = np.subtract(b1, bar_to_subtract[2])
                        bar_subtracted = np.subtract(b2, bar_to_subtract[3])

                        bars_subtracted[i] = bar_subtracted

                    print "residual evap"
                    print bars_subtracted

                    bars_subtracted = bars_subtracted / variable_data_5

                    print variable_data_1
                    print variable_data_2
                    print variable_data_3
                    print bars_subtracted
                    print variable_data_4
                    print variable_data_5
                    #print color

                    p5 = plt.bar(number_of_bars, bars_subtracted, bar_width, bottom=np.array(variable_data_1)+np.array(variable_data_2)+np.array(variable_data_3), color=color)

                    for bar in p5:
                        bar.set_edgecolor("black")
                        bar.set_linewidth(0.5)
                        bar.set_linestyle("-")

        if include_dynamic_prescribed == 'yes':
            model_strings_for_plot = np.append(model_strings_for_plot, "Dynamic Vegetation")
            model_strings_for_plot = np.append(model_strings_for_plot, "Prescribed Vegetation")

        if ensemble == 'yes':
            model_strings_for_plot = np.append(model_strings_for_plot, "Ensemble")
        if number_of_reanalysis > 0:
            model_strings_for_plot = np.append(model_strings_for_plot, reanalysis_strings_for_plot)

        count = 0

        for i in model_strings_for_plot:

            if i in ['bcc-csm1-1-m', 'BNU-ESM', 'GFDL-HIRAM-C360', 'IPSL-CM5A-MR', 'Dynamic Vegetation']:

                model_string = str(i) + "*"

                model_strings_for_plot[count] = model_string

                count+=1

                #model_strings_for_plot[count] = r"$\bf{"+str(i)+"}$"

            if i in ['CanAM4', 'CNRM-CM5', 'GISS-E2-R', 'inmcm4', 'MIROC5', 'MRI-AGCM3-2S', 'NorESM1-M', 'Prescribed Vegetation', 'Ensemble', 'GLEAM']:

                #model_strings_for_plot[count] = r"$\it{"+str(i)+"}$"

                count+=1

        plt.ylabel(bar_season+" "+bar_y_axis_title, fontsize=9)
        plt.xticks(number_of_bars, model_strings_for_plot, rotation=45, ha='right', fontsize=9)

        plt.ylim(lower_y_lim, upper_y_lim)

        for i in list_of_variables:
            if i == 'tran':
                list_of_variables = [i.replace('tran', 'Transpiration') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'evspsblsoi':
                list_of_variables = [i.replace('evspsblsoi', 'Bare Soil Evaporation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'evspsblveg':
                list_of_variables = [i.replace('evspsblveg', 'Canopy Evaporation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'evaporation':
                list_of_variables = [i.replace('evaporation', 'Evaporation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'nrad':
                list_of_variables = [i.replace('nrad', 'Surface Net Downward Radiation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'hfls':
                list_of_variables = [i.replace('hfls', 'Evaporation') for i in list_of_variables]
        for i in list_of_variables:
            if i == 'hfss':
                list_of_variables = [i.replace('hfss', 'Surface Upward Sensible Heat Flux') for i in list_of_variables]

        if len(variables_to_subtract) == 0:

            plt.legend((p4[0], p3[0], p2[0], p1[0]), list_of_variables[::-1], fontsize=9)

        if len(variables_to_subtract) > 0:

            if divide_four_variables_by_fifth != 'yes':

                list_of_variables.insert(3, string)

                print list_of_variables
                print list_of_variables[::-1]
                plt.legend((p4[0], p5[0], p3[0], p2[0], p1[0]), list_of_variables[::-1], fontsize=8)

            if divide_four_variables_by_fifth == 'yes':

                list_of_variables.insert(3, string)

                print list_of_variables

                list_of_variables = list_of_variables[:-1]

                print list_of_variables[::-1]

                plt.legend((p4[0], p5[0], p3[0], p2[0], p1[0]), list_of_variables[::-1], fontsize=8)


        plt.gcf().subplots_adjust(bottom=0.3)

        if bar_season == 'DJF':
            bar_season = 'S1'
        if bar_season == 'MAM':
            bar_season = 'S2'
        if bar_season == 'JJA':
            bar_season = 'S3'
        if bar_season == 'SON':
            bar_season = 'S4'

        print "Saving figure"
        fig.savefig("Bar_Graph_"+bar_season, bbox_inches='tight')


def plot_seasonal_cycle(seasonal_cycle_models_array, seasonal_cycle_ensemble_array, seasonal_cycle_ensemble_std_array, seasonal_cycle_ensemble_top_array, seasonal_cycle_ensemble_bottom_array, seasonal_cycle_dv_array, seasonal_cycle_pv_array, seasonal_cycle_reanalysis_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim, y_tick_interval, number_of_variables, list_of_variables, variables_to_add, fill_between_lines, variables_to_subtract, legend_in_plot, number_of_models, list_of_models):

    """Before looping through all the models, set up the figure to plot to."""

#-------------------------------------------------------------------------------------

    if subplot_multiple_variables != 'yes':
        print "Setting up figure"
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        if duplicate_seasonal_cycle != 'yes':
            objects = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan')
            ax1.set_xlim([0, 12])
        if duplicate_seasonal_cycle == 'yes':
            objects = ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J')
            ax1.set_xlim([0, 24])
        x_pos = np.arange(len(objects))
        plt.xticks(x_pos, objects)
        ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
        ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')
        ax1.patch.set_visible(False)

    #-------------------------------------------------------------------------------------

        print number_of_variables
        if number_of_variables == 1:
            print "Figure set up"

            """Load the data from the model file paths into a cube. Constrain the input years"""
            """Print the model ID, length of time dimension, and first and last model dates."""

            """Plot model seasonal cycles."""

            if len(model_strings_for_plot) == 0:
                pass
            else:

                for i in np.arange(0, len(seasonal_cycle_models_array)):

                    if duplicate_seasonal_cycle != 'yes':

                        seasonal_cycle_model = seasonal_cycle_models_array[i]

                    if duplicate_seasonal_cycle == 'yes':
                        seasonal_cycle_model = seasonal_cycle_models_array[i]
                        seasonal_cycle_model_2 = seasonal_cycle_models_array[i]
                        seasonal_cycle_model_2 = seasonal_cycle_model_2[1:]
                        seasonal_cycle_model = np.append(seasonal_cycle_model, seasonal_cycle_model_2)


                    line_colour = model_line_colours[i]
                    model_id = model_strings_for_plot[i]
                    print model_id

                    print x_pos
                    print seasonal_cycle_model

                    ax1.plot(x_pos, seasonal_cycle_model, zorder=1, linestyle='-', color=line_colour, label = str(model_id))
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("-")
                    if include_legend == 'yes':
                        legend = plt.legend(handles, labels, loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9, handlelength=2.5)

                """Plot model ensemble mean seasonal cycle."""

            if len(ensemble_string_for_plot) == 0:
                pass

            else:

                if duplicate_seasonal_cycle != 'yes':

                    ensemble_mean_seasonal_cycle_models = seasonal_cycle_ensemble_array

                if duplicate_seasonal_cycle == 'yes':
                    ensemble_mean_seasonal_cycle_models = seasonal_cycle_ensemble_array
                    ensemble_mean_seasonal_cycle_models_2 = seasonal_cycle_ensemble_array
                    ensemble_mean_seasonal_cycle_models_2 = seasonal_cycle_ensemble_array[1:]
                    ensemble_mean_seasonal_cycle_models = np.append(ensemble_mean_seasonal_cycle_models, ensemble_mean_seasonal_cycle_models_2)

                ensemble_string = ensemble_string_for_plot[0]
                print ensemble_string

                ax1.plot(x_pos, ensemble_mean_seasonal_cycle_models, zorder=2, linestyle='-', linewidth=4.0, color='black', label = str(ensemble_string))
                handles, labels = ax1.get_legend_handles_labels()
                handles[-1].set_linestyle("-")
                if include_legend == 'yes':
                    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

            """Plot reanalysis seasonal cycles."""

            if len(reanalysis_strings_for_plot) == 0:
                pass
            else:

                for i in np.arange(0, len(seasonal_cycle_reanalysis_array)):

                    if duplicate_seasonal_cycle != 'yes':

                        seasonal_cycle_reanalysis = seasonal_cycle_reanalysis_array[i]

                    if duplicate_seasonal_cycle == 'yes':
                        seasonal_cycle_reanalysis = seasonal_cycle_reanalysis_array[i]
                        seasonal_cycle_reanalysis_2 = seasonal_cycle_reanalysis_array[i]
                        seasonal_cycle_reanalysis_2 = seasonal_cycle_reanalysis_2[1:]
                        seasonal_cycle_reanalysis = np.append(seasonal_cycle_reanalysis, seasonal_cycle_reanalysis_2)

                    line_colour = reanalysis_line_colours[i]
                    reanalysis_id = reanalysis_strings_for_plot[i]
                    print reanalysis_id

                    ax1.plot(x_pos, seasonal_cycle_reanalysis, zorder=1, linestyle='--', color = line_colour, label = str(reanalysis_id))
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("--")
                    if include_legend == 'yes':
                        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

            if include_dynamic_prescribed == 'yes':

                if duplicate_seasonal_cycle != 'yes':

                    seasonal_cycle_dv = seasonal_cycle_dv_array

                if duplicate_seasonal_cycle == 'yes':
                    print "hi5"
                    print seasonal_cycle_dv_array

                    seasonal_cycle_dv = seasonal_cycle_dv_array
                    seasonal_cycle_dv_2 = seasonal_cycle_dv_array
                    seasonal_cycle_dv_2 = seasonal_cycle_dv_2[1:]
                    seasonal_cycle_dv = np.append(seasonal_cycle_dv, seasonal_cycle_dv_2)

                if duplicate_seasonal_cycle != 'yes':

                    seasonal_cycle_pv = seasonal_cycle_pv_array

                if duplicate_seasonal_cycle == 'yes':
                    seasonal_cycle_pv = seasonal_cycle_pv_array
                    seasonal_cycle_pv_2 = seasonal_cycle_pv_array
                    seasonal_cycle_pv_2 = seasonal_cycle_pv_2[1:]
                    seasonal_cycle_pv = np.append(seasonal_cycle_pv, seasonal_cycle_pv_2)

                print "hi6"
                line = ax1.plot(x_pos, seasonal_cycle_dv, linestyle='-', linewidth=4.0, color='forestgreen', label = "Best Models")
                handles, labels = ax1.get_legend_handles_labels()
                handles[-1].set_linestyle("-")
                if include_legend == 'yes':
                    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

                line = ax1.plot(x_pos, seasonal_cycle_pv, linestyle='-', linewidth=4.0, color='saddlebrown', label = "Worst Models")
                handles, labels = ax1.get_legend_handles_labels()
                handles[-1].set_linestyle("-")
                if include_legend == 'yes':
                    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

            print lower_y_lim
            print upper_y_lim+y_tick_interval
            print y_tick_interval
            print np.arange(lower_y_lim, upper_y_lim+y_tick_interval, y_tick_interval)
            print "hi3"
            print lower_y_lim
            print upper_y_lim
            plt.ylim(lower_y_lim, upper_y_lim)
            plt.yticks(np.arange(lower_y_lim, upper_y_lim+y_tick_interval, y_tick_interval))

            if variable == 'pr':
                plt.ylabel('Precipitation (mm $\mathregular{day^{-1}}$)')
            if variable == 'hfls':
                plt.ylabel('Surface Upward Latent Heat Flux (W $\mathregular{m^{-2}}$)')
            if variable == 'evaporation':
                plt.ylabel('Evaporation (mm $\mathregular{day^{-1}}$)')
            if variable == 'hfss':
                plt.ylabel('Surface Upward Sensible Heat Flux (W $\mathregular{m^{-2}}$)', fontsize=10)
            if variable == 'evap_fraction':
                plt.ylabel('Evaporative Fraction')
            if variable == 'nrad':
                plt.ylabel('Surface Net Downward Radiation (W $\mathregular{m^{-2}}$)')
            if variable == 'mrsos':
                plt.ylabel('Soil Moisture Content of Upper Layer (mm)')
            if variable == 'mrso':
                plt.ylabel('Soil Moisture Content (mm)')
            if variable == 'tran':
                if unit_plot == 'mm day-1':
                    plt.ylabel('Transpiration (mm $\mathregular{day^{-1}}$)')
                if unit_plot == 'W m-2':
                    plt.ylabel('Transpiration (W $\mathregular{m^{-2}}$)')
            if variable == 'evspsbl':
                if unit_plot == 'mm day-1':
                    plt.ylabel('Evaporation (mm $\mathregular{day^{-1}}$)')
                if unit_plot == 'W m-2':
                    plt.ylabel('Evaporation (W $\mathregular{m^{-2}}$)')
            if variable == 'evspsblsoi':
                if unit_plot == 'mm day-1':
                    plt.ylabel('Bare Soil Evaporation (mm $\mathregular{day^{-1}}$)')
                if unit_plot == 'W m-2':
                        plt.ylabel('Bare Soil Evaporation (W $\mathregular{m^{-2}}$)')
            if variable == 'evspsblveg':
                if unit_plot == 'mm day-1':
                    plt.ylabel('Canopy Evaporation (mm $\mathregular{day^{-1}}$)')
                if unit_plot == 'W m-2':
                    plt.ylabel('Canopy Evaporation (W $\mathregular{m^{-2}}$)')
            if variable == 'prveg':
                plt.ylabel('Precipitation Intercepted by Canopy (mm $\mathregular{day^{-1}}$)')
            if variable == 'mrros':
                plt.ylabel('Surface Runoff Flux (mm $\mathregular{day^{-1}}$)')
            if variable == 'lai':
                plt.ylabel('Leaf Area Index')
            if variable == 'mrro':
                plt.ylabel('Total Runoff (mm $\mathregular{day^{-1}}$)')
            if variable == 'vpd':
                plt.ylabel('Vapour Pressure Deficit (kPa)')
            if variable == 'hurs':
                plt.ylabel('Surface Relative Humidity (%)')
            if variable == 'tas':
                plt.ylabel('Surface Air Temperature ($^\circ$C)')
            if variable == 'sa':
                plt.ylabel('Soil Moisture Accumulation (mm $\mathregular{day^{-1}}$)')
            if variable == 'swc_anom':
                plt.ylabel('Soil Water Content Anomaly (mm)')
            if variable == 'cancap':
                plt.ylabel('Canopy Storage Capacity (mm)')

            """Save the figure."""
            print "Saving figure"
            if include_legend == 'yes':
                fig.savefig("Seasonal_Cycle_"+variable+".png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)
            if include_legend != 'yes':
                fig.savefig("Seasonal_Cycle_"+variable+".png", bbox_inches='tight', dpi=600)
            print "Plot done."

        #-------------------------------------------------------------------------------------

        if number_of_variables > 1:

            #JUST REANALYSIS AND ENSEMBLE FOR NOW

            print "hi8"
            print "Figure set up"
            print len(variables_to_add)
            # for i in variables_to_add:
            #     for j in np.arange(0, len(i)):
            #         print i[j]

            """Load the data from the model file paths into a cube. Constrain the input years"""
            """Print the model ID, length of time dimension, and first and last model dates."""

            """Plot model seasonal cycles."""

            #-------------------------------------------------------------------------------------

            if len(reanalysis_strings_for_plot) > 0:

                variable_number = 0
                line_number = 0
                line_count = 0

                fill_between_lines_data = seasonal_cycle_to_add = np.zeros((2, 13))

                #-------------------------------------------------------------------------------------

                for i in np.arange(0, len(seasonal_cycle_reanalysis_array)):
                    print "hi9"
                    seasonal_cycle_reanalysis_variable_data = seasonal_cycle_reanalysis_array[variable_number]

                    name_of_variable = list_of_variables[variable_number]

                    #['saddlebrown', 'dodgerblue', 'forestgreen']

                    if name_of_variable == 'evaporation':
                        reanalysis_string = 'Evaporation'
                        color = 'black'
                        zorder = 2

                    if name_of_variable == 'pr':
                        reanalysis_string = 'Precipitation'
                        color = 'blue'
                        zorder = 3

                    if name_of_variable == 'mrros':
                        reanalysis_string = 'Surface Runoff'
                        color = 'green'
                        zorder = 1

                    if name_of_variable == 'mrro':
                        reanalysis_string = 'Runoff'
                        color = 'green'
                        zorder = 1

                    if name_of_variable == 'tran':
                        reanalysis_string = 'Transpiration'
                        color = 'dodgerblue'
                        zorder = 1

                    if name_of_variable == 'evspsblsoi':
                        reanalysis_string = 'Bare Soil Evaporation'
                        color = 'saddlebrown'
                        zorder = 1

                    if name_of_variable == 'evspsblveg':
                        reanalysis_string = 'Canopy Evaporation'
                        color = 'forestgreen'
                        zorder = 1

                    line = ax1.plot(x_pos, seasonal_cycle_reanalysis_variable_data, zorder=zorder, linestyle='-', linewidth=2.0, color=color, label = reanalysis_string)
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("-")

                    if legend_in_plot != 'yes':
                        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                    if legend_in_plot == 'yes':
                        legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                    variable_number +=1
                    line_number +=1

                    if len(fill_between_lines) > 0:

                        if line_number in fill_between_lines:
                            fill_between_lines_data[line_count] = seasonal_cycle_reanalysis_variable_data
                            line_count +=1

                        if line_number not in fill_between_lines:
                            line[0].remove()

                #-------------------------------------------------------------------------------------

                if len(variables_to_add) > 0:

                    indices_to_slice_seasonal_cycle = np.zeros((len(variables_to_add), 2))
                    variables_to_add_count = 0
                    for i in variables_to_add:
                        print i

                        if i == ['pr', 'evaporation']:
                            reanalysis_string = 'AMIP Precipitation + Evaporation'
                            color = 'green'
                        if i == ['evaporation', 'mrros']:
                            reanalysis_string = 'Evaporation + Surface Runoff'
                            color = 'red'
                        if i == ['evaporation', 'mrro']:
                            reanalysis_string = 'Evaporation + Runoff'
                            color = 'red'

                        indices = []

                        for j in i:
                            indices = np.append(indices, list_of_variables.index(j))
                        print indices
                        indices_to_slice_seasonal_cycle[variables_to_add_count] = indices
                    print indices_to_slice_seasonal_cycle

                    for i in indices_to_slice_seasonal_cycle:
                        print i
                        seasonal_cycle_to_add = np.zeros((len(i), 13))
                        print seasonal_cycle_to_add.shape
                        count = 0
                        for j in i:
                            seasonal_cycle_to_add_data = seasonal_cycle_reanalysis_array[int(j)]
                            seasonal_cycle_to_add[count] = seasonal_cycle_to_add_data
                            count +=1
                        print seasonal_cycle_to_add

                        seasonal_cycle_added = seasonal_cycle_to_add.sum(axis = 0)

                        line = ax1.plot(x_pos, seasonal_cycle_added, zorder=2, linestyle='-', linewidth=2.0, color=color, label = reanalysis_string)
                        handles, labels = ax1.get_legend_handles_labels()
                        handles[-1].set_linestyle("-")

                        if legend_in_plot != 'yes':
                            legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                        if legend_in_plot == 'yes':
                            legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                        if len(fill_between_lines) > 0:

                            if line_number in fill_between_lines:
                                fill_between_lines_data[line_count] = seasonal_cycle_added
                                line_count +=1

                            if line_number not in fill_between_lines:
                                line[0].remove()

                #-------------------------------------------------------------------------------------

                if len(variables_to_subtract) > 0:
                    print variables_to_subtract

                    variables_to_subtract_count = 0
                    for i in variables_to_subtract:
                        print i

                        if i == ['pr', 'evspsblsoi', 'evspsblveg', 'mrro']:
                            reanalysis_string = 'Available Water'
                            color = 'blue'

                        if i == ['pr', 'mrro']:
                            reanalysis_string = 'Available Water'
                            color = 'blue'

                        if i == ['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']:
                            reanalysis_string = 'Residual Evaporation'
                            color = 'red'

                        if i == ['pr', 'tran', 'evspsblveg', 'mrro']:
                            reanalysis_string = 'Available Water'
                            color = 'blue'

                        indices_to_slice_seasonal_cycle = np.zeros((1, len(i)))

                        indices = []
                        for j in i:
                            print j
                            print list_of_variables.index(j)
                            indices = np.append(indices, list_of_variables.index(j))

                        print indices_to_slice_seasonal_cycle.shape
                        print variables_to_subtract_count

                        indices_to_slice_seasonal_cycle[variables_to_subtract_count] = indices

                    print indices_to_slice_seasonal_cycle

                    for i in indices_to_slice_seasonal_cycle:
                        print i
                        seasonal_cycle_to_subtract = np.zeros((len(i), 13))
                        print seasonal_cycle_to_subtract.shape
                        count = 0
                        for j in i:
                            seasonal_cycle_to_subtract_data = seasonal_cycle_reanalysis_array[int(j)]
                            seasonal_cycle_to_subtract[count] = seasonal_cycle_to_subtract_data
                            count +=1
                        print seasonal_cycle_to_subtract

                        if len(i) == 4:

                            s1 = np.subtract(seasonal_cycle_to_subtract[0], seasonal_cycle_to_subtract[1])
                            s2 = np.subtract(s1, seasonal_cycle_to_subtract[2])
                            seasonal_cycle_subtracted = np.subtract(s2, seasonal_cycle_to_subtract[3])

                            print seasonal_cycle_subtracted

                        line = ax1.plot(x_pos, seasonal_cycle_subtracted, zorder=2, linestyle='-', linewidth=2.0, color=color, label = reanalysis_string)
                        handles, labels = ax1.get_legend_handles_labels()
                        handles[-1].set_linestyle("-")

                        if legend_in_plot != 'yes':
                            legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                        if legend_in_plot == 'yes':
                            legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                        line_number +=1

                        if len(fill_between_lines) > 0:

                            if line_number in fill_between_lines:
                                fill_between_lines_data[line_count] = seasonal_cycle_subtracted
                                line_count +=1

                            if line_number not in fill_between_lines:
                                line[0].remove()

                #-------------------------------------------------------------------------------------

                if len(fill_between_lines) > 0:

                    print fill_between_lines_data

                    y1 = fill_between_lines_data[0]
                    y2 = fill_between_lines_data[1]

                    ax1.fill_between(x_pos, y1, y2, where=y2 <= y1, facecolor='lightseagreen', interpolate=True, label = 'Soil Moisture Accumulation')
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("-")

                    if legend_in_plot != 'yes':
                        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                    if legend_in_plot == 'yes':
                        legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                    ax1.fill_between(x_pos, y1, y2, where=y2 >= y1, facecolor='saddlebrown', interpolate=True, label = 'Soil Moisture Depletion')
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("-")

                    if legend_in_plot != 'yes':
                        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                    if legend_in_plot == 'yes':
                        legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

            #-------------------------------------------------------------------------------------

            if len(ensemble_string_for_plot) > 0:

                variable_number = 0
                line_number = 0
                line_count = 0

                fill_between_lines_data = seasonal_cycle_to_add = np.zeros((2, 13))

                for i in np.arange(0, len(seasonal_cycle_ensemble_array)):
                    print "hi9"
                    seasonal_cycle_ensemble_variable_data = seasonal_cycle_ensemble_array[variable_number]

                    name_of_variable = list_of_variables[variable_number]

                    #['saddlebrown', 'dodgerblue', 'forestgreen']

                    if name_of_variable == 'evaporation':
                        ensemble_string = 'Evaporation'
                        color = 'red'
                        zorder = 2

                    if name_of_variable == 'pr':
                        ensemble_string = 'Precipitation'
                        color = 'blue'
                        zorder = 3

                    if name_of_variable == 'mrros':
                        ensemble_string = 'Surface Runoff'
                        color = 'green'
                        zorder = 1

                    if name_of_variable == 'mrro':
                        ensemble_string = 'Runoff'
                        color = 'green'
                        zorder = 1

                    if name_of_variable == 'tran':
                        ensemble_string = 'Transpiration'
                        color = 'dodgerblue'
                        zorder = 1

                    if name_of_variable == 'evspsblsoi':
                        ensemble_string = 'Bare Soil Evaporation'
                        color = 'saddlebrown'
                        zorder = 1

                    if name_of_variable == 'evspsblveg':
                        ensemble_string = 'Canopy Evaporation'
                        color = 'forestgreen'
                        zorder = 1

                    line = ax1.plot(x_pos, seasonal_cycle_ensemble_variable_data, zorder=zorder, linestyle='-', linewidth=2.0, color=color, label = ensemble_string)
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("-")

                    if legend_in_plot != 'yes':
                        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                    if legend_in_plot == 'yes':
                        legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                    if len(fill_between_lines) > 0:

                        if line_number in fill_between_lines:
                            fill_between_lines_data[line_count] = seasonal_cycle_ensemble_variable_data
                            line_count +=1

                        if line_number not in fill_between_lines:
                            line[0].remove()

                    variable_number +=1
                    line_number +=1

                if len(variables_to_add) > 0:

                    indices_to_slice_seasonal_cycle = np.zeros((len(variables_to_add), 2))
                    variables_to_add_count = 0
                    for i in variables_to_add:
                        print i

                        if i == ['pr', 'evaporation']:
                            ensemble_string = 'AMIP Precipitation + Evaporation'
                            color = 'green'
                        if i == ['evaporation', 'mrros']:
                            ensemble_string = 'Evaporation + Surface Runoff'
                            color = 'red'
                        if i == ['evaporation', 'mrro']:
                            ensemble_string = 'Evaporation + Runoff'
                            color = 'red'

                        indices = []

                        for j in i:
                            indices = np.append(indices, list_of_variables.index(j))
                        print indices
                        indices_to_slice_seasonal_cycle[variables_to_add_count] = indices
                    print indices_to_slice_seasonal_cycle

                    for i in indices_to_slice_seasonal_cycle:
                        print i
                        seasonal_cycle_to_add = np.zeros((len(i), 13))
                        print seasonal_cycle_to_add.shape
                        count = 0
                        for j in i:
                            seasonal_cycle_to_add_data = seasonal_cycle_ensemble_array[int(j)]
                            seasonal_cycle_to_add[count] = seasonal_cycle_to_add_data
                            count +=1
                        print seasonal_cycle_to_add

                        seasonal_cycle_added = seasonal_cycle_to_add.sum(axis = 0)

                        line = ax1.plot(x_pos, seasonal_cycle_added, zorder=2, linestyle='-', linewidth=2.0, color=color, label = ensemble_string)
                        handles, labels = ax1.get_legend_handles_labels()
                        handles[-1].set_linestyle("-")

                        if legend_in_plot != 'yes':
                            legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                        if legend_in_plot == 'yes':
                            legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                        if len(fill_between_lines) > 0:

                            if line_number in fill_between_lines:
                                fill_between_lines_data[line_count] = seasonal_cycle_added
                                line_count +=1

                            if line_number not in fill_between_lines:
                                line[0].remove()

                if len(variables_to_subtract) > 0:
                    print variables_to_subtract


                    variables_to_subtract_count = 0
                    for i in variables_to_subtract:
                        print i

                        if i == ['pr', 'evspsblsoi', 'evspsblveg', 'mrro']:
                            ensemble_string = 'Available Water'
                            color = 'blue'

                        if i == ['pr', 'mrro']:
                            ensemble_string = 'Available Water'
                            color = 'blue'

                        if i == ['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']:
                            ensemble_string = 'Residual Evaporation'
                            color = 'red'

                        if i == ['pr', 'tran', 'evspsblveg', 'mrro']:
                            ensemble_string = 'Available Water'
                            color = 'blue'

                        indices_to_slice_seasonal_cycle = np.zeros((1, len(i)))

                        indices = []
                        for j in i:
                            print j
                            print list_of_variables.index(j)
                            indices = np.append(indices, list_of_variables.index(j))

                        print indices_to_slice_seasonal_cycle.shape
                        print variables_to_subtract_count

                        indices_to_slice_seasonal_cycle[variables_to_subtract_count] = indices

                    print indices_to_slice_seasonal_cycle

                    for i in indices_to_slice_seasonal_cycle:
                        print i
                        seasonal_cycle_to_subtract = np.zeros((len(i), 13))
                        print seasonal_cycle_to_subtract.shape
                        count = 0
                        for j in i:
                            seasonal_cycle_to_subtract_data = seasonal_cycle_ensemble_array[int(j)]
                            seasonal_cycle_to_subtract[count] = seasonal_cycle_to_subtract_data
                            count +=1
                        print seasonal_cycle_to_subtract

                        if len(i) == 4:

                            s1 = np.subtract(seasonal_cycle_to_subtract[0], seasonal_cycle_to_subtract[1])
                            s2 = np.subtract(s1, seasonal_cycle_to_subtract[2])
                            seasonal_cycle_subtracted = np.subtract(s2, seasonal_cycle_to_subtract[3])

                            print seasonal_cycle_subtracted

                        if len(i) == 2:

                            print "hi18"
                            print seasonal_cycle_to_subtract
                            seasonal_cycle_subtracted = np.subtract(seasonal_cycle_to_subtract[0], seasonal_cycle_to_subtract[1])
                            print seasonal_cycle_subtracted

                        line = ax1.plot(x_pos, seasonal_cycle_subtracted, zorder=2, linestyle='-', linewidth=2.0, color=color, label = ensemble_string)
                        handles, labels = ax1.get_legend_handles_labels()
                        handles[-1].set_linestyle("-")

                        if legend_in_plot != 'yes':
                            legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                        if legend_in_plot == 'yes':
                            legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                        line_number +=1

                        print line_number

                        if len(fill_between_lines) > 0:

                            if line_number in fill_between_lines:
                                fill_between_lines_data[line_count] = seasonal_cycle_subtracted
                                line_count +=1

                            if line_number not in fill_between_lines:
                                line[0].remove()

                if len(fill_between_lines) > 0:
                    print "hi19"
                    y1 = fill_between_lines_data[0]
                    y2 = fill_between_lines_data[1]

                    print y1
                    print y2

                    for i in variables_to_subtract:
                        if i == ['pr', 'evspsblsoi', 'evspsblveg', 'mrro']:
                            y1 = fill_between_lines_data[1]
                            y2 = fill_between_lines_data[0]

                        if i == ['pr', 'mrro']:
                            y1 = fill_between_lines_data[1]
                            y2 = fill_between_lines_data[0]

                    print y1
                    print y2

                    ax1.fill_between(x_pos, y1, y2, where=y2 <= y1, facecolor='lightseagreen', interpolate=True, label = 'Soil Moisture Accumulation')
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("-")

                    if legend_in_plot != 'yes':
                        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                    if legend_in_plot == 'yes':
                        legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                    ax1.fill_between(x_pos, y1, y2, where=y2 >= y1, facecolor='saddlebrown', interpolate=True, label = 'Soil Moisture Depletion')
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("-")

                    if legend_in_plot != 'yes':
                        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)
                    if legend_in_plot == 'yes':
                        legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)

                for i in variables_to_subtract:
                    if i == ['pr', 'evspsblsoi', 'evspsblveg', 'mrro']:
                        handles, labels = ax1.get_legend_handles_labels()
                        handles[0], handles[1] = handles[1], handles[0]
                        labels[0], labels[1] = labels[1], labels[0]
                        legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)
                        print handles
                        print labels

                    if i == ['pr', 'mrro']:
                        handles, labels = ax1.get_legend_handles_labels()
                        handles[0], handles[1] = handles[1], handles[0]
                        labels[0], labels[1] = labels[1], labels[0]
                        legend = plt.legend(handles, labels, fontsize=9, handlelength=2.5)
                        print handles
                        print labels

            print lower_y_lim
            print upper_y_lim+y_tick_interval
            print y_tick_interval
            print np.arange(lower_y_lim, upper_y_lim+y_tick_interval, y_tick_interval)
            plt.ylim(lower_y_lim, upper_y_lim)
            plt.yticks(np.arange(lower_y_lim, upper_y_lim+y_tick_interval, y_tick_interval))

            if list_of_variables[0] == 'pr' or list_of_variables[0] == 'evaporation':
                plt.ylabel('Evaporation (mm $\mathregular{day^{-1}}$)')

            """Save the figure."""
            print "Saving figure"
            fig.savefig("Seasonal_Cycle_"+variable+"_north.png", bbox_inches='tight', dpi=600)
            print "Plot done."

    #-------------------------------------------------------------------------------------

    if subplot_multiple_variables == 'yes' and number_of_variables >= 1:

        if subplot_columns == 4 and subplot_rows == 2:
            fig = plt.figure(figsize=(8,6))

        if subplot_columns == 4 and subplot_rows == 3:
            fig = plt.figure(figsize=(10,8))

        if subplot_columns == 4 and subplot_rows == 4:
            fig = plt.figure(figsize=(10, 10))

        if subplot_columns == 4 and subplot_rows == 5:
            fig = plt.figure(figsize=(10, 12))

        if subplot_columns == 1 and subplot_rows == 2:
            fig = plt.figure(figsize=(10, 12))

        #def plot_seasonal_cycle(seasonal_cycle_models_array, seasonal_cycle_ensemble_array, seasonal_cycle_reanalysis_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim, y_tick_interval, number_of_variables, list_of_variables, variables_to_add, fill_between_lines, variables_to_subtract, legend_in_plot):

        variable_number = 0

        for i in range(len(list_of_variables)):

            ax = plt.subplot(subplot_rows, subplot_columns, variable_number+1)
            ax.set_xlim([0, 12])
            ax.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
            ax.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')
            ax.patch.set_visible(False)

            objects = ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J')
            x_pos = np.arange(len(objects))
            plt.xticks(x_pos, objects, fontsize=7)

            if list_of_variables[i] == 'vpd':
                lower_y_lim = 0.0
                upper_y_lim = 3.2
                y_tick_interval = 0.4
                plt.ylabel('Vapour Pressure Deficit (kPa)', fontsize=7)
                line_colour = 'steelblue'

            if list_of_variables[i] == 'tas':
                lower_y_lim = 10
                upper_y_lim = 35
                y_tick_interval = 5
                plt.ylabel('Air Temperature ($^\circ$C )', fontsize=7)
                line_colour = 'darkorange'

            if list_of_variables[i] == 'hurs':
                lower_y_lim = 10
                upper_y_lim = 100
                y_tick_interval = 10
                plt.ylabel('Relative Humidity (%)', fontsize=7)
                line_colour = 'deepskyblue'

            if list_of_variables[i] == 'tran':
                lower_y_lim = 0.0
                upper_y_lim = 5.5
                y_tick_interval = 0.5
                plt.ylabel('Transpiration (mm $\mathregular{day^{-1}}$)', fontsize=7)
                line_colour = 'dodgerblue'

            if list_of_variables[i] == 'evspsblveg':
                lower_y_lim = 0.0
                upper_y_lim = 5.5
                y_tick_interval = 0.5
                plt.ylabel('Canopy Evaporation (mm $\mathregular{day^{-1}}$)', fontsize=7)
                line_colour = 'forestgreen'

            if list_of_variables[i] == 'evspsblsoi':
                lower_y_lim = 0.0
                upper_y_lim = 5.5
                y_tick_interval = 0.5
                plt.ylabel('Bare Soil Evaporation (mm $\mathregular{day^{-1}}$)', fontsize=7)
                line_colour = 'saddlebrown'

            if list_of_variables[i] == 'evaporation':
                lower_y_lim = 0.0
                upper_y_lim = 5.5
                y_tick_interval = 0.5
                plt.ylabel('Evaporation (mm $\mathregular{day^{-1}}$)', fontsize=7)
                line_colour = 'black'

            if list_of_variables[i] == 'pr':
                lower_y_lim = 0.0
                upper_y_lim = 13.0
                y_tick_interval = 1.0
                plt.ylabel('Precipitation (mm $\mathregular{day^{-1}}$)', fontsize=7)
                line_colour = 'blue'

            if list_of_variables[i] == 'lai':
                lower_y_lim = 0
                upper_y_lim = 8
                y_tick_interval = 1
                plt.ylabel('Leaf Area Index', fontsize=7)
                line_colour = 'darkgreen'

            if list_of_variables[i] == 'nrad':
                lower_y_lim = 80
                upper_y_lim = 190
                y_tick_interval = 10
                plt.ylabel('Net Downward Radiation (W $\mathregular{m^{-2}}$)', fontsize=7)
                line_colour = 'gold'

            if list_of_variables[i] == 'mrsos':
                lower_y_lim = 0
                upper_y_lim = 45
                y_tick_interval = 5
                plt.ylabel('Upper Layer Soil Moisture (mm)', fontsize=7)
                line_colour = 'blueviolet'

            if list_of_variables[i] == 'ws':
                lower_y_lim = 0.0
                upper_y_lim = 6.4
                y_tick_interval = 0.8
                plt.ylabel('Wind Speed (m s $\mathregular{^{-1}}$)', fontsize=7)
                line_colour = 'peru'

            if list_of_variables[i] == 'sa':
                lower_y_lim = -4.0
                upper_y_lim = 4.8
                y_tick_interval = 0.8
                plt.ylabel('Soil Moisture Accumulation (mm $\mathregular{day^{-1}}$)', fontsize=7)
                line_colour = 'darkcyan'

            if list_of_variables[i] == 'hfls':
                lower_y_lim = 0.0
                upper_y_lim = 160.0
                y_tick_interval = 20.0
                plt.ylabel('Upward Latent Heat Flux (W $\mathregular{m^{-2}}$)', fontsize=7)
                line_colour = 'navy'

            if list_of_variables[i] == 'hfss':
                lower_y_lim = 0.0
                upper_y_lim = 110.0
                y_tick_interval = 10.0
                plt.ylabel('Upward Sensible Heat Flux (W $\mathregular{m^{-2}}$)', fontsize=7)
                line_colour = 'darkred'

            if list_of_variables[i] == 'mrro':
                lower_y_lim = -1.0
                upper_y_lim = 9.0
                y_tick_interval = 1.0
                plt.ylabel('Total Runoff (mm $\mathregular{day^{-1}}$)', fontsize=7)
                line_colour = 'purple'

            if list_of_variables[i] == 'cancap':
                lower_y_lim = 0
                upper_y_lim = 2
                y_tick_interval = 0.2
                plt.ylabel('Canopy Storage Capacity (mm)', fontsize=7)
                line_colour = 'olive'

            plt.ylim(lower_y_lim, upper_y_lim)
            plt.yticks(np.arange(lower_y_lim, upper_y_lim+y_tick_interval, y_tick_interval), fontsize=7)

            if models == 'yes':

                for j in range(len(list_of_models)):

                    line_colour = model_line_colours[j]

                    seasonal_cycle_model_variable_data = seasonal_cycle_models_array[i,:,j]

                    model_id = model_strings_for_plot[j]

                    line = ax.plot(x_pos, seasonal_cycle_model_variable_data, zorder=1, linestyle='-', color=line_colour, label = str(model_id))

                    handles, labels = ax.get_legend_handles_labels()
                    handles[-1].set_linestyle("-")
                    legend = plt.legend(handles, labels, loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9, handlelength=2.5)

            if ensemble == 'yes':

                seasonal_cycle_ensemble_variable_data = seasonal_cycle_ensemble_array[i]

                #seasonal_cycle_ensemble_variable_std_data = seasonal_cycle_ensemble_std_array[i]

                seasonal_cycle_ensemble_variable_top_data = seasonal_cycle_ensemble_top_array[i]

                seasonal_cycle_ensemble_variable_bottom_data = seasonal_cycle_ensemble_bottom_array[i]

                print list_of_variables[i]
                print seasonal_cycle_ensemble_variable_data
                print "max"
                print seasonal_cycle_ensemble_array[i] + seasonal_cycle_ensemble_top_array[i]
                print "min"
                print seasonal_cycle_ensemble_array[i] - seasonal_cycle_ensemble_bottom_array[i]

                #line = ax.plot(x_pos, seasonal_cycle_ensemble_variable_data, linestyle='-', linewidth=1.5, color=line_colour, zorder=2)
                line = ax.plot(x_pos, seasonal_cycle_ensemble_variable_data, zorder=2, linestyle='-', linewidth=4.0, color='black', label = str(ensemble_string))
                handles, labels = ax.get_legend_handles_labels()
                handles[-1].set_linestyle("-")
                legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

                """
                line2 = ax.errorbar(x_pos, seasonal_cycle_ensemble_variable_data, yerr=(seasonal_cycle_ensemble_variable_bottom_data, seasonal_cycle_ensemble_variable_top_data), fmt='none', ecolor='black', elinewidth=0.5, capsize=2, capthick=1.25, linestyle = '--', zorder=1)

                line2[-1][0].set_linestyle('--')

                xticks, xticklabels = plt.xticks()

                xmin = (2*xticks[0] - xticks[1])

                xmax = (2*xticks[-1] - xticks[-2])

                ax.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='off', top='off')

                plt.xlim(xmin, xmax)
                plt.xticks(xticks)
                """

            if reanalysis == 'yes':

                seasonal_cycle_reanalysis_variable_data = seasonal_cycle_reanalysis_array[i]

                line = ax.plot(x_pos, seasonal_cycle_reanalysis_variable_data, linestyle='--', linewidth=1.5, color=line_colour, zorder=3)

            if include_dynamic_prescribed == 'yes':
                print "hi6"
                seasonal_cycle_dv_variable_data = seasonal_cycle_dv_array[i]
                line = ax.plot(x_pos, seasonal_cycle_dv_variable_data, linestyle='-', linewidth=4.0, color='forestgreen')

                seasonal_cycle_pv_variable_data = seasonal_cycle_pv_array[i]
                line = ax.plot(x_pos, seasonal_cycle_pv_variable_data, linestyle='-', linewidth=4.0, color='saddlebrown')

            variable_number +=1

        fig.tight_layout()

        if subplot_columns == 4 and subplot_rows == 2:

            #plt.subplots_adjust(left=0.05, right=0.95, wspace=0.8)
            plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

        if subplot_columns == 4 and subplot_rows == 3:

            plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

        if subplot_columns == 4 and subplot_rows == 4:

            plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

        if subplot_columns == 4 and subplot_rows == 5:

            plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

        if subplot_columns == 1 and subplot_rows == 2:

            plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

        print "saving final figure"
        fig.savefig("test_multiple_models.png", bbox_inches='tight', dpi=600)
        plt.close()
        print "plot done"


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "CanAM4/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["CNRM-CM5/"]
    list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    list_of_reanalysis = ["modis"]

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['lai']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'no'
    # list_dynamic_veg_models = ['bcc-csm1-1/', 'bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C180/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-LR/', 'IPSL-CM5A-MR/', 'IPSL-CM5B-LR/']
    # list_prescribed_veg_models = ['CanAM4/', 'CNRM-CM5/', 'GISS-E2-R/', 'inmcm4/', 'MIROC5/', "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/']
    #list_prescribed_veg_models = ['CanAM4/', 'CNRM-CM5/', 'GISS-E2-R/', 'inmcm4/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/"]
    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ['BNU-ESM/', 'MIROC5/', "MRI-AGCM3-2S/"]

    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

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
    lower_y_lim = 0.0
    upper_y_lim = 10.0
    y_tick_interval = 1.0
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0

    cmap = 'rainbow'
    models = 'no'
    ensemble = 'no'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'


    # ------------------------------------------------------------------------------------------------------------------------------------------

    time_range = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)

    # ------------------------------------------------------------------------------------------------------------------------------------------

    """Set up arrays for bar chart."""

    print len(list_of_variables)
    print len(list_of_models)

    if ensemble == 'no':

        barchart_array_djf = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis))))
        barchart_array_mam = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis))))
        barchart_array_jja = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis))))
        barchart_array_son = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis))))
        barchart_array_mar = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis))))
        barchart_array_jul = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis))))
        barchart_array_nov = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis))))

    if ensemble == 'yes':

        barchart_array_djf = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_mam = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_jja = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_son = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_mar = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_jul = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_nov = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))


    if ensemble == 'yes' and include_dynamic_prescribed == 'yes':

        barchart_array_djf = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+3)))
        barchart_array_mam = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+3)))
        barchart_array_jja = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+3)))
        barchart_array_son = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+3)))
        barchart_array_mar = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+3)))
        barchart_array_jul = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+3)))
        barchart_array_nov = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+3)))


    """For each variable, calculate the model, ensemble and reanalysis seasonal cycle arrays. """

    """Define first variable number."""

    variable_number = 0

    seasonal_cycle_models_multiple_variables_array = np.zeros((len(list_of_variables), 13, len(list_of_models)))
    seasonal_cycle_ensemble_multiple_variables_array = np.zeros((len(list_of_variables), 13))
    seasonal_cycle_ensemble_multiple_variables_std_array = np.zeros((len(list_of_variables), 13))
    seasonal_cycle_ensemble_multiple_variables_top_array = np.zeros((len(list_of_variables), 13))
    seasonal_cycle_ensemble_multiple_variables_bottom_array = np.zeros((len(list_of_variables), 13))
    seasonal_cycle_pv_multiple_variables_array = np.zeros((len(list_of_variables), 13))
    seasonal_cycle_dv_multiple_variables_array = np.zeros((len(list_of_variables), 13))
    seasonal_cycle_reanalysis_multiple_variables_array = np.zeros((len(list_of_variables), 13))

    list_of_models2 = list_of_models

    for i in list_of_variables:

        seasonal_cycle_models_array = []
        seasonal_cycle_ensemble_array = []
        seasonal_cycle_ensemble_std_array = []
        seasonal_cycle_ensemble_top_array = []
        seasonal_cycle_ensemble_bottom_array = []
        seasonal_cycle_dv_array = []
        seasonal_cycle_pv_array = []
        seasonal_cycle_reanalysis_array = []
        model_strings_for_plot = []
        ensemble_string_for_plot = []
        reanalysis_strings_for_plot = []

        variable = i
        print variable

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

                        """Multiprocess calculations for each month."""
                        manager = multiprocessing.Manager()
                        array = manager.dict()
                        jobs = []
                        for i in np.arange(1, 13, 1):
                            p = multiprocessing.Process(target=seasonal_cycle_calculations, args=(i, array))
                            jobs.append(p)
                            p.start()
                        for process in jobs:
                            process.join()
                        model_seasonal_cycle = array.values()
                        print model_id
                        print model_seasonal_cycle
                        model_strings_for_plot = np.append(model_strings_for_plot, model_id)
                        jan_value = model_seasonal_cycle[0]
                        model_seasonal_cycle.append(jan_value)
                        seasonal_cycle_models_array.append(model_seasonal_cycle)

    # ------------------------------------------------------------------------------------------------------------------------------------------

            """ENSEMBLE MEAN."""

            if ensemble == 'yes':

                if variable == 'lai':

                    #list_to_remove_lai = ['BNU-ESM/', 'CanAM4/', 'CNRM-CM5/', 'GISS-E2-R/', 'inmcm4/', 'MRI-AGCM3-2H/', 'MRI-AGCM3-2S/', 'MRI-CGCM3/']

                    list_to_remove_lai = []

                    list_to_remove_lai = set(list_to_remove_lai)

                    list_models_lai = [x for x in list_of_models if x not in list_to_remove_lai]

                    print list_models_lai

                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_models_lai, model_type, variable)

                    if include_dynamic_prescribed == 'yes':

                        list_models_lai_dv = [x for x in list_models_lai if x in list_dynamic_veg_models]

                        model_file_paths_ensemble_dv = model_file_paths_ensemble_func(list_models_lai_dv, model_type, variable)

                        list_models_lai_pv = [x for x in list_models_lai if x in list_prescribed_veg_models]

                        model_file_paths_ensemble_pv = model_file_paths_ensemble_func(list_models_lai_pv, model_type, variable)

                        print model_file_paths_ensemble
                        print model_file_paths_ensemble_pv
                        print model_file_paths_ensemble_dv

                if variable == 'cancap':

                    list_to_remove_cancap = ['BNU-ESM/', 'CanAM4/', 'CNRM-CM5/', 'GISS-E2-R/', 'inmcm4/', 'MRI-AGCM3-2H/', 'MRI-AGCM3-2S/', 'MRI-CGCM3/']

                    list_to_remove_cancap = set(list_to_remove_cancap)

                    list_models_cancap = [x for x in list_of_models if x not in list_to_remove_cancap]

                    print list_models_cancap

                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_models_cancap, model_type, variable)

                    if include_dynamic_prescribed == 'yes':

                        list_models_cancap_dv = [x for x in list_models_cancap if x in list_dynamic_veg_models]

                        model_file_paths_ensemble_dv = model_file_paths_ensemble_func(list_models_cancap_dv, model_type, variable)

                        list_models_cancap_pv = [x for x in list_models_cancap if x in list_prescribed_veg_models]

                        model_file_paths_ensemble_pv = model_file_paths_ensemble_func(list_models_cancap_pv, model_type, variable)

                        print model_file_paths_ensemble
                        print model_file_paths_ensemble_pv
                        print model_file_paths_ensemble_dv

                if variable == 'mrsos':

                    list_to_remove_mrsos = ['IPSL-CM5A-LR/', 'IPSL-CM5A-MR/', 'IPSL-CM5B-LR/']

                    list_to_remove_mrsos = set(list_to_remove_mrsos)

                    list_models_mrsos = [x for x in list_of_models if x not in list_to_remove_mrsos]

                    print list_models_mrsos

                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_models_mrsos, model_type, variable)

                    if include_dynamic_prescribed == 'yes':

                        list_models_mrsos_dv = [x for x in list_models_mrsos if x in list_dynamic_veg_models]

                        model_file_paths_ensemble_dv = model_file_paths_ensemble_func(list_models_mrsos_dv, model_type, variable)

                        list_models_mrsos_pv = [x for x in list_models_mrsos if x in list_prescribed_veg_models]

                        model_file_paths_ensemble_pv = model_file_paths_ensemble_func(list_models_mrsos_pv, model_type, variable)

                        print model_file_paths_ensemble
                        print model_file_paths_ensemble_pv
                        print model_file_paths_ensemble_dv

                if variable == 'vpd':

                    list_to_remove_vpd = ['GFDL-HIRAM-C180/', 'GFDL-HIRAM-C360/', 'GISS-E2-R/', 'IPSL-CM5A-LR/', 'IPSL-CM5A-MR/', 'IPSL-CM5B-LR/', 'MIROC5/', 'MRI-AGCM3-2H/', 'MRI-AGCM3-2S/', 'MRI-CGCM3/']

                    list_to_remove_vpd = set(list_to_remove_vpd)

                    list_models_vpd = [x for x in list_of_models if x not in list_to_remove_vpd]

                    print list_models_vpd

                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_models_vpd, model_type, variable)

                    if include_dynamic_prescribed == 'yes':

                        list_models_vpd_dv = [x for x in list_models_vpd if x in list_dynamic_veg_models]

                        model_file_paths_ensemble_dv = model_file_paths_ensemble_func(list_models_vpd_dv, model_type, variable)

                        list_models_vpd_pv = [x for x in list_models_vpd if x in list_prescribed_veg_models]

                        model_file_paths_ensemble_pv = model_file_paths_ensemble_func(list_models_vpd_pv, model_type, variable)

                        print model_file_paths_ensemble
                        print model_file_paths_ensemble_pv
                        print model_file_paths_ensemble_dv

                if variable == 'hurs':

                    list_to_remove_hurs = ['GFDL-HIRAM-C180/', 'GFDL-HIRAM-C360/', 'GISS-E2-R/', 'IPSL-CM5A-LR/', 'IPSL-CM5A-MR/', 'IPSL-CM5B-LR/', 'MIROC5/', 'MRI-AGCM3-2H/', 'MRI-AGCM3-2S/', 'MRI-CGCM3/']

                    list_to_remove_hurs = set(list_to_remove_hurs)

                    list_models_hurs = [x for x in list_of_models if x not in list_to_remove_hurs]

                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_models_hurs, model_type, variable)

                    if include_dynamic_prescribed == 'yes':

                        list_models_hurs_dv = [x for x in list_models_hurs if x in list_dynamic_veg_models]

                        model_file_paths_ensemble_dv = model_file_paths_ensemble_func(list_models_hurs_dv, model_type, variable)

                        list_models_hurs_pv = [x for x in list_models_hurs if x in list_prescribed_veg_models]

                        model_file_paths_ensemble_pv = model_file_paths_ensemble_func(list_models_hurs_pv, model_type, variable)

                        print model_file_paths_ensemble
                        print model_file_paths_ensemble_pv
                        print model_file_paths_ensemble_dv

                if variable not in ['mrsos', 'lai', 'cancap', 'hurs', 'vpd']:

                    """Extract the regridded model file paths for the ensemble mean."""
                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

                    model_file_paths_ensemble_dv = model_file_paths_ensemble_func(list_dynamic_veg_models, model_type, variable)

                    model_file_paths_ensemble_pv = model_file_paths_ensemble_func(list_prescribed_veg_models, model_type, variable)

                    print model_file_paths_ensemble
                    print model_file_paths_ensemble_dv
                    print model_file_paths_ensemble_pv

                #---------------------------------------------------------------------------------------------------

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

                seasonal_cycle_array_for_ensemble = []

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
                        for i in np.arange(1, 13, 1):
                            p = multiprocessing.Process(target=seasonal_cycle_calculations, args=(i, array))
                            jobs.append(p)
                            p.start()
                        for process in jobs:
                            process.join()
                        cube_seasonal_cycle_for_ensemble = array.values()
                        print model_id
                        print cube_seasonal_cycle_for_ensemble
                        jan_value = cube_seasonal_cycle_for_ensemble[0]
                        cube_seasonal_cycle_for_ensemble.append(jan_value)
                        seasonal_cycle_array_for_ensemble.append(cube_seasonal_cycle_for_ensemble)

                seasonal_cycle_ensemble_array = np.mean(seasonal_cycle_array_for_ensemble, axis = 0)

                seasonal_cycle_ensemble_std_array = np.std(seasonal_cycle_array_for_ensemble, axis = 0)

                seasonal_cycle_ensemble_max_array = np.amax(seasonal_cycle_array_for_ensemble, axis = 0)

                seasonal_cycle_ensemble_min_array = np.amin(seasonal_cycle_array_for_ensemble, axis = 0)

                seasonal_cycle_ensemble_top_array = seasonal_cycle_ensemble_max_array - seasonal_cycle_ensemble_array

                seasonal_cycle_ensemble_bottom_array = seasonal_cycle_ensemble_array - seasonal_cycle_ensemble_min_array

                print "ensemble"
                print seasonal_cycle_ensemble_array

                ensemble_string_for_plot = np.append(ensemble_string_for_plot, "Ensemble")

                #---------------------------------------------------------------------------------------------------

                if include_dynamic_prescribed == 'yes':

                #---------------------------------------------------------------------------------------------------

                    """DYNAMIC VEG ENSEMBLE MEAN."""

                    model_file_paths_ensemble = model_file_paths_ensemble_dv

                    print model_file_paths_ensemble
                    print "hi7"

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

                    seasonal_cycle_array_for_ensemble_dv = []

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
                            for i in np.arange(1, 13, 1):
                                p = multiprocessing.Process(target=seasonal_cycle_calculations, args=(i, array))
                                jobs.append(p)
                                p.start()
                            for process in jobs:
                                process.join()
                            cube_seasonal_cycle_for_ensemble = array.values()
                            print model_id
                            print cube_seasonal_cycle_for_ensemble
                            jan_value = cube_seasonal_cycle_for_ensemble[0]
                            cube_seasonal_cycle_for_ensemble.append(jan_value)
                            seasonal_cycle_array_for_ensemble_dv.append(cube_seasonal_cycle_for_ensemble)

                    seasonal_cycle_dv_array = np.mean(seasonal_cycle_array_for_ensemble_dv, axis = 0)
                    print "dynamic veg"
                    print seasonal_cycle_dv_array

                    #---------------------------------------------------------------------------------------------------

                    """PRESCRIBED VEG ENSEMBLE MEAN."""

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

                    seasonal_cycle_array_for_ensemble_pv = []

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
                            for i in np.arange(1, 13, 1):
                                p = multiprocessing.Process(target=seasonal_cycle_calculations, args=(i, array))
                                jobs.append(p)
                                p.start()
                            for process in jobs:
                                process.join()
                            cube_seasonal_cycle_for_ensemble = array.values()
                            print model_id
                            print cube_seasonal_cycle_for_ensemble
                            jan_value = cube_seasonal_cycle_for_ensemble[0]
                            cube_seasonal_cycle_for_ensemble.append(jan_value)
                            seasonal_cycle_array_for_ensemble_pv.append(cube_seasonal_cycle_for_ensemble)

                    seasonal_cycle_pv_array = np.mean(seasonal_cycle_array_for_ensemble_pv, axis = 0)
                    print "prescribed veg"
                    print seasonal_cycle_pv_array


        # ------------------------------------------------------------------------------------------------------------------------------------------

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
                    for i in np.arange(1, 13, 1):
                        p = multiprocessing.Process(target=seasonal_cycle_calculations, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cube_seasonal_cycle_reanalysis = array.values()
                    print reanalysis_id
                    print cube_seasonal_cycle_reanalysis
                    jan_value = cube_seasonal_cycle_reanalysis[0]
                    cube_seasonal_cycle_reanalysis.append(jan_value)
                    reanalysis_strings_for_plot = np.append(reanalysis_strings_for_plot, reanalysis_id)
                    seasonal_cycle_reanalysis_array.append(cube_seasonal_cycle_reanalysis)

        #------------------------------------------------------------------------------------------------------------------------------------------

        """test to see if arrays are correct size."""

        # if len(list_of_models) > 0:
        #     print np.shape(seasonal_cycle_models_array)
        # if ensemble == 'yes':
        #     print np.shape(seasonal_cycle_ensemble_array)
        # if len(list_of_reanalysis) > 0:
        #     print np.shape(seasonal_cycle_reanalysis_array)

        print('')

        #------------------------------------------------------------------------------------------------------------------------------------------

        """PLOT SEASONAL CYCLE FOR EACH VARIABLE IF SELECTED."""

        if len(list_of_models) > 0:
            print model_strings_for_plot
        if ensemble == "yes":
            print ensemble_string_for_plot
        if len(list_of_reanalysis) > 0:
            print reanalysis_strings_for_plot

        print('')

        number_of_models = len(list_of_models)
        number_of_reanalysis = len(list_of_reanalysis)
        number_of_variables = len(list_of_variables)

        print "hi32"
        print seasonal_cycle_models_array
        print seasonal_cycle_ensemble_array
        print seasonal_cycle_reanalysis_array

        if plot == 'yes' and subplot_multiple_variables != 'yes':

            if number_of_variables == 1:

                model_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[0]
                reanalysis_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[1]
                plot_seasonal_cycle(seasonal_cycle_models_array, seasonal_cycle_ensemble_array, seasonal_cycle_ensemble_std_array, seasonal_cycle_ensemble_top_array, seasonal_cycle_ensemble_bottom_array, seasonal_cycle_dv_array, seasonal_cycle_pv_array, seasonal_cycle_reanalysis_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim, y_tick_interval, number_of_variables, list_of_variables, variables_to_add, fill_between_lines, variables_to_subtract, legend_in_plot, number_of_models, list_of_models)
                print('')


    #------------------------------------------------------------------------------------------------------------------------------------------

        """PRINT VALUES FOR SEASONAL CYCLE."""

        djf_values = []
        mam_values = []
        jja_values = []
        son_values = []
        mar_values = []
        jul_values = []
        nov_values = []
        list_all_datasets = []

        if len(list_of_models) > 0:
            list_all_datasets = np.append(list_all_datasets, list_of_models)

        if include_dynamic_prescribed == 'yes':
            list_all_datasets = np.append(list_all_datasets, "Dynamic Veg")
            list_all_datasets = np.append(list_all_datasets, "Prescribed Veg")

        if ensemble == "yes":
            list_all_datasets = np.append(list_all_datasets, "Ensemble")
        if len(list_of_reanalysis) > 0:
            list_all_datasets = np.append(list_all_datasets, list_of_reanalysis)

        print "hi3"
        print model_strings_for_plot

        if len(list_of_models) == 0:
            model_number = -1

        if len(list_of_models) > 0:

            if models == 'no':
                pass

            else:

                for i in np.arange(0, len(list_of_models)):

                    model_number = i

                    model_id = model_strings_for_plot[i]
                    print model_id
                    model_seasonal_cycle = seasonal_cycle_models_array[i]

                    model_seasonal_cycle = [float('%.3f' % i) for i in model_seasonal_cycle]

                    print model_seasonal_cycle

                    djf_value = (model_seasonal_cycle[11] + model_seasonal_cycle[0] + model_seasonal_cycle[1])/float(3.0)
                    mam_value = (model_seasonal_cycle[2] + model_seasonal_cycle[3] + model_seasonal_cycle[4])/float(3.0)
                    jja_value = (model_seasonal_cycle[5] + model_seasonal_cycle[6] + model_seasonal_cycle[7])/float(3.0)
                    son_value = (model_seasonal_cycle[8] + model_seasonal_cycle[9] + model_seasonal_cycle[10])/float(3.0)
                    mar_value = model_seasonal_cycle[2]
                    jul_value = model_seasonal_cycle[6]
                    nov_value = model_seasonal_cycle[10]

                    print "%.3f, %.3f, %.3f, %.3f" % (djf_value, mam_value, jja_value, son_value)

                    print('')

                    djf_values = np.append(djf_values, djf_value)
                    mam_values = np.append(mam_values, mam_value)
                    jja_values = np.append(jja_values, jja_value)
                    son_values = np.append(son_values, son_value)
                    mar_values = np.append(mar_values, mar_value)
                    jul_values = np.append(jul_values, jul_value)
                    nov_values = np.append(nov_values, nov_value)

                    barchart_array_djf[variable_number, model_number] = djf_value
                    barchart_array_mam[variable_number, model_number] = mam_value
                    barchart_array_jja[variable_number, model_number] = jja_value
                    barchart_array_son[variable_number, model_number] = son_value
                    barchart_array_mar[variable_number, model_number] = mar_value
                    barchart_array_jul[variable_number, model_number] = jul_value

                    seasonal_cycle_models_multiple_variables_array[variable_number,:,model_number] = model_seasonal_cycle

                    print variable_number
                    print model_number

        if ensemble == 'yes':

            if include_dynamic_prescribed == 'yes' and models == 'yes':

                dv_seasonal_cycle =  [float('%.3f' % i) for i in seasonal_cycle_dv_array]

                seasonal_cycle_dv_multiple_variables_array[variable_number] = dv_seasonal_cycle

                pv_seasonal_cycle =  [float('%.3f' % i) for i in seasonal_cycle_pv_array]

                seasonal_cycle_pv_multiple_variables_array[variable_number] = pv_seasonal_cycle


                djf_value_dv = (dv_seasonal_cycle[11] + dv_seasonal_cycle[0] + dv_seasonal_cycle[1])/float(3.0)
                mam_value_dv = (dv_seasonal_cycle[2] + dv_seasonal_cycle[3] + dv_seasonal_cycle[4])/float(3.0)
                jja_value_dv = (dv_seasonal_cycle[5] + dv_seasonal_cycle[6] + dv_seasonal_cycle[7])/float(3.0)
                son_value_dv = (dv_seasonal_cycle[8] + dv_seasonal_cycle[9] + dv_seasonal_cycle[10])/float(3.0)
                mar_value_dv = dv_seasonal_cycle[2]
                jul_value_dv = dv_seasonal_cycle[6]
                nov_value_dv = dv_seasonal_cycle[10]

                djf_value_pv = (pv_seasonal_cycle[11] + pv_seasonal_cycle[0] + pv_seasonal_cycle[1])/float(3.0)
                mam_value_pv = (pv_seasonal_cycle[2] + pv_seasonal_cycle[3] + pv_seasonal_cycle[4])/float(3.0)
                jja_value_pv = (pv_seasonal_cycle[5] + pv_seasonal_cycle[6] + pv_seasonal_cycle[7])/float(3.0)
                son_value_pv = (pv_seasonal_cycle[8] + pv_seasonal_cycle[9] + pv_seasonal_cycle[10])/float(3.0)
                mar_value_pv = pv_seasonal_cycle[2]
                jul_value_pv = pv_seasonal_cycle[6]
                nov_value_pv = pv_seasonal_cycle[10]

                djf_values = np.append(djf_values, djf_value_dv)
                mam_values = np.append(mam_values, mam_value_dv)
                jja_values = np.append(jja_values, jja_value_dv)
                son_values = np.append(son_values, son_value_dv)
                mar_values = np.append(mar_values, mar_value_dv)
                jul_values = np.append(jul_values, jul_value_dv)
                nov_values = np.append(nov_values, nov_value_dv)

                djf_values = np.append(djf_values, djf_value_pv)
                mam_values = np.append(mam_values, mam_value_pv)
                jja_values = np.append(jja_values, jja_value_pv)
                son_values = np.append(son_values, son_value_pv)
                mar_values = np.append(mar_values, mar_value_pv)
                jul_values = np.append(jul_values, jul_value_pv)
                nov_values = np.append(nov_values, nov_value_pv)

                model_number = model_number+1

                barchart_array_djf[variable_number, model_number] = djf_value_dv
                barchart_array_mam[variable_number, model_number] = mam_value_dv
                barchart_array_jja[variable_number, model_number] = jja_value_dv
                barchart_array_son[variable_number, model_number] = son_value_dv
                barchart_array_mar[variable_number, model_number] = mar_value_dv
                barchart_array_jul[variable_number, model_number] = jul_value_dv
                barchart_array_nov[variable_number, model_number] = nov_value_dv

                model_number = model_number+1

                barchart_array_djf[variable_number, model_number] = djf_value_pv
                barchart_array_mam[variable_number, model_number] = mam_value_pv
                barchart_array_jja[variable_number, model_number] = jja_value_pv
                barchart_array_son[variable_number, model_number] = son_value_pv
                barchart_array_mar[variable_number, model_number] = mar_value_pv
                barchart_array_jul[variable_number, model_number] = jul_value_pv
                barchart_array_nov[variable_number, model_number] = nov_value_pv

                print "Best Models"
                print dv_seasonal_cycle
                print "Worst Models"
                print pv_seasonal_cycle

            print "Ensemble"
            ensemble_seasonal_cycle = [float('%.3f' % i) for i in seasonal_cycle_ensemble_array]
            ensemble_std_seasonal_cycle = [float('%.3f' % i) for i in seasonal_cycle_ensemble_std_array]
            ensemble_top_seasonal_cycle = [float('%.3f' % i) for i in seasonal_cycle_ensemble_top_array]
            ensemble_bottom_seasonal_cycle = [float('%.3f' % i) for i in seasonal_cycle_ensemble_bottom_array]
            print variable
            print ensemble_seasonal_cycle

            djf_value = (ensemble_seasonal_cycle[11] + ensemble_seasonal_cycle[0] + ensemble_seasonal_cycle[1])/float(3.0)
            mam_value = (ensemble_seasonal_cycle[2] + ensemble_seasonal_cycle[3] + ensemble_seasonal_cycle[4])/float(3.0)
            jja_value = (ensemble_seasonal_cycle[5] + ensemble_seasonal_cycle[6] + ensemble_seasonal_cycle[7])/float(3.0)
            son_value = (ensemble_seasonal_cycle[8] + ensemble_seasonal_cycle[9] + ensemble_seasonal_cycle[10])/float(3.0)
            mar_value = ensemble_seasonal_cycle[2]
            jul_value = ensemble_seasonal_cycle[6]
            nov_value = ensemble_seasonal_cycle[10]

            print "%.3f, %.3f, %.3f, %.3f" % (djf_value, mam_value, jja_value, son_value)

            #print('')

            djf_values = np.append(djf_values, djf_value)
            mam_values = np.append(mam_values, mam_value)
            jja_values = np.append(jja_values, jja_value)
            son_values = np.append(son_values, son_value)
            mar_values = np.append(mar_values, mar_value)
            jul_values = np.append(jul_values, jul_value)
            nov_values = np.append(nov_values, nov_value)

            if models == 'yes':
                model_number = model_number+1
            if models == 'no':
                model_number = 0

            barchart_array_djf[variable_number, model_number] = djf_value
            barchart_array_mam[variable_number, model_number] = mam_value
            barchart_array_jja[variable_number, model_number] = jja_value
            barchart_array_son[variable_number, model_number] = son_value
            barchart_array_mar[variable_number, model_number] = mar_value
            barchart_array_jul[variable_number, model_number] = jul_value

            print variable_number

            seasonal_cycle_ensemble_multiple_variables_array[variable_number] = ensemble_seasonal_cycle
            seasonal_cycle_ensemble_multiple_variables_std_array[variable_number] = ensemble_std_seasonal_cycle
            seasonal_cycle_ensemble_multiple_variables_top_array[variable_number] = ensemble_top_seasonal_cycle
            seasonal_cycle_ensemble_multiple_variables_bottom_array[variable_number] = ensemble_bottom_seasonal_cycle


        if len(list_of_reanalysis) > 0:

            for i in np.arange(0, len(seasonal_cycle_reanalysis_array)):

                reanalysis_number = i

                reanalysis_id = reanalysis_strings_for_plot[i]
                print reanalysis_id
                reanalysis_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                reanalysis_seasonal_cycle = [float('%.3f' % i) for i in reanalysis_seasonal_cycle]
                print reanalysis_seasonal_cycle

                djf_value = (reanalysis_seasonal_cycle[11] + reanalysis_seasonal_cycle[0] + reanalysis_seasonal_cycle[1])/float(3.0)
                mam_value = (reanalysis_seasonal_cycle[2] + reanalysis_seasonal_cycle[3] + reanalysis_seasonal_cycle[4])/float(3.0)
                jja_value = (reanalysis_seasonal_cycle[5] + reanalysis_seasonal_cycle[6] + reanalysis_seasonal_cycle[7])/float(3.0)
                son_value = (reanalysis_seasonal_cycle[8] + reanalysis_seasonal_cycle[9] + reanalysis_seasonal_cycle[10])/float(3.0)
                mar_value = reanalysis_seasonal_cycle[2]
                jul_value = reanalysis_seasonal_cycle[6]
                nov_value = reanalysis_seasonal_cycle[10]

                print "%.3f, %.3f, %.3f, %.3f" % (djf_value, mam_value, jja_value, son_value)

                print('')

                djf_values = np.append(djf_values, djf_value)
                mam_values = np.append(mam_values, mam_value)
                jja_values = np.append(jja_values, jja_value)
                son_values = np.append(son_values, son_value)
                mar_values = np.append(mar_values, mar_value)
                jul_values = np.append(jul_values, jul_value)
                nov_values = np.append(nov_values, nov_value)

                print mar_values

                reanalysis_number = reanalysis_number + 1 + model_number

                barchart_array_djf[variable_number, reanalysis_number] = djf_value
                barchart_array_mam[variable_number, reanalysis_number] = mam_value
                barchart_array_jja[variable_number, reanalysis_number] = jja_value
                barchart_array_son[variable_number, reanalysis_number] = son_value
                barchart_array_mar[variable_number, reanalysis_number] = mar_value
                barchart_array_jul[variable_number,  reanalysis_number] = jul_value

                print barchart_array_mar

                seasonal_cycle_reanalysis_multiple_variables_array[variable_number] = reanalysis_seasonal_cycle

        list_all_datasets_input = [str(i) for i in list_all_datasets]

        print list_all_datasets_input

        print "DJF VALUES"

        print djf_values

        djf_values, list_all_datasets = (list(i) for i in zip(*sorted(zip(djf_values, list_all_datasets_input))))
        djf_values = [float('%.3f' % i) for i in djf_values]

        print djf_values[::-1]
        print list_all_datasets[::-1]
        print list_all_datasets

        print('')

        print "MAM VALUES"

        print mam_values
        mam_values, list_all_datasets = (list(i) for i in zip(*sorted(zip(mam_values, list_all_datasets_input))))
        mam_values = [float('%.3f' % i) for i in mam_values]
        print mam_values[::-1]
        print list_all_datasets[::-1]
        print list_all_datasets

        print('')

        print "JJA VALUES"

        print jja_values

        jja_values, list_all_datasets = (list(i) for i in zip(*sorted(zip(jja_values, list_all_datasets_input))))
        jja_values = [float('%.3f' % i) for i in jja_values]
        print jja_values[::-1]
        print list_all_datasets[::-1]
        print list_all_datasets

        print('')

        print "SON VALUES"

        print son_values
        son_values, list_all_datasets = (list(i) for i in zip(*sorted(zip(son_values, list_all_datasets_input))))
        son_values = [float('%.3f' % i) for i in son_values]
        print son_values[::-1]
        print list_all_datasets[::-1]
        print list_all_datasets

        print('')

        print "Mar VALUES"

        print mar_values
        # print list_all_datasets_input
        # mar_values, list_all_datasets = (list(i) for i in zip(*sorted(zip(mar_values, list_all_datasets_input))))
        # mar_values = [float('%.3f' % i) for i in mar_values]
        # print mar_values[::-1]
        # print list_all_datasets[::-1]
        # print list_all_datasets

        print('')

        print "Jul VALUES"

        print jul_values
        jul_values, list_all_datasets = (list(i) for i in zip(*sorted(zip(jul_values, list_all_datasets_input))))
        jul_values = [float('%.3f' % i) for i in jul_values]
        print jul_values[::-1]
        print list_all_datasets[::-1]
        print list_all_datasets

        print "Nov VALUES"

        print nov_values
        # nov_values, list_all_datasets = (list(i) for i in zip(*sorted(zip(nov_values, list_all_datasets_input))))
        # nov_values = [float('%.3f' % i) for i in nov_values]
        # print nov_values[::-1]
        # print list_all_datasets[::-1]
        # print list_all_datasets

        print "Nov - Mar VALUES"

        nov_mar_values = np.array(nov_values) - np.array(mar_values)
        print nov_mar_values
        print list_all_datasets_input


        print('')


        """
        print mam_values
        print jja_values
        print son_values
        print list_all_datasets
        """

        if rmse == 'yes':

            months = [1,2,3,4,9,10,11,12]
            months[:] = [i - 1 for i in months]

            evap_rmse = []
            pr_rmse = []
            list_of_rmse_models_evap = []
            list_of_rmse_models_precip = []

            if len(list_of_models) > 0:

                for i in np.arange(0, len(list_of_models)):

                    model_number = i
                    model_id = model_strings_for_plot[i]
                    model_seasonal_cycle = seasonal_cycle_models_array[i]
                    model_seasonal_cycle = [float('%.3f' % i) for i in model_seasonal_cycle]

                    if len(list_of_reanalysis) > 0:

                        array_storing_reanalysis = np.zeros((2, 12))

                        for i in np.arange(0, len(list_of_reanalysis)):

                            reanalysis_number = i

                            reanalysis_id = reanalysis_strings_for_plot[i]

                            if variable == 'evaporation' and reanalysis_id == 'LandFlux-EVAL':

                                landfluxeval_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                                landfluxeval_seasonal_cycle = [float('%.3f' % i) for i in landfluxeval_seasonal_cycle]

                                model_seasonal_cycle = model_seasonal_cycle[:-1]
                                landfluxeval_seasonal_cycle = landfluxeval_seasonal_cycle[:-1]

                                model_seasonal_cycle = np.asarray(model_seasonal_cycle)
                                landfluxeval_seasonal_cycle = np.asarray(landfluxeval_seasonal_cycle)

                                print model_id
                                print model_seasonal_cycle[months]
                                print landfluxeval_seasonal_cycle[months]

                                rmse_val = np.sqrt(((np.array(model_seasonal_cycle[months]) - np.array(landfluxeval_seasonal_cycle[months])) ** 2).mean())
                                print rmse_val

                                evap_rmse = np.append(evap_rmse, rmse_val)
                                list_of_rmse_models_evap = np.append(list_of_rmse_models_evap, model_id)

                            if variable == 'pr':

                                if reanalysis_id == 'CHIRPS':
                                    chirps_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                                    chirps_seasonal_cycle = [float('%.3f' % i) for i in chirps_seasonal_cycle]
                                    chirps_seasonal_cycle = chirps_seasonal_cycle[:-1]
                                    array_storing_reanalysis[0] = chirps_seasonal_cycle

                                if reanalysis_id == 'GPCC':
                                    gpcc_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                                    gpcc_seasonal_cycle = [float('%.3f' % i) for i in gpcc_seasonal_cycle]
                                    gpcc_seasonal_cycle = gpcc_seasonal_cycle[:-1]
                                    array_storing_reanalysis[1] = gpcc_seasonal_cycle

                                pr_reference_array = np.mean(array_storing_reanalysis, axis = 0)

                                pr_ref_seasonal_cycle = [float('%.3f' % i) for i in pr_reference_array]

                        if variable == 'pr':
                            model_seasonal_cycle = model_seasonal_cycle[:-1]
                            model_seasonal_cycle = np.asarray(model_seasonal_cycle)
                            pr_ref_seasonal_cycle = np.asarray(pr_ref_seasonal_cycle)
                            print model_id
                            print model_seasonal_cycle
                            print pr_ref_seasonal_cycle
                            rmse_val = np.sqrt(((np.array(model_seasonal_cycle[months]) - np.array(pr_ref_seasonal_cycle[months])) ** 2).mean())
                            print rmse_val
                            pr_rmse = np.append(pr_rmse, rmse_val)
                            list_of_rmse_models_precip = np.append(list_of_rmse_models_precip, model_id)

            if ensemble == 'yes':

                ensemble_seasonal_cycle = [float('%.3f' % i) for i in seasonal_cycle_ensemble_array]
                ensemble_seasonal_cycle = ensemble_seasonal_cycle[:-1]

                if len(list_of_reanalysis) > 0:

                    array_storing_reanalysis = np.zeros((2, 12))

                    for i in np.arange(0, len(seasonal_cycle_reanalysis_array)):

                        reanalysis_number = i

                        reanalysis_id = reanalysis_strings_for_plot[i]

                        if variable == 'evaporation' and reanalysis_id == 'LandFlux-EVAL':

                            landfluxeval_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                            landfluxeval_seasonal_cycle = [float('%.3f' % i) for i in landfluxeval_seasonal_cycle]
                            landfluxeval_seasonal_cycle = landfluxeval_seasonal_cycle[:-1]

                            print "ensemble"

                            ensemble_seasonal_cycle = np.asarray(ensemble_seasonal_cycle)
                            landfluxeval_seasonal_cycle = np.asarray(landfluxeval_seasonal_cycle)

                            print ensemble_seasonal_cycle[months]
                            print landfluxeval_seasonal_cycle[months]

                            rmse_val = np.sqrt(((np.array(ensemble_seasonal_cycle[months]) - np.array(landfluxeval_seasonal_cycle[months])) ** 2).mean())

                            evap_rmse = np.append(evap_rmse, rmse_val)
                            list_of_rmse_models_evap = np.append(list_of_rmse_models_evap, "Ensemble")

                            if variable == 'pr':

                                if reanalysis_id == 'CHIRPS':
                                    chirps_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                                    chirps_seasonal_cycle = [float('%.3f' % i) for i in chirps_seasonal_cycle]
                                    chirps_seasonal_cycle = chirps_seasonal_cycle[:-1]
                                    array_storing_reanalysis[0] = chirps_seasonal_cycle

                                if reanalysis_id == 'GPCC':
                                    gpcc_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                                    gpcc_seasonal_cycle = [float('%.3f' % i) for i in gpcc_seasonal_cycle]
                                    gpcc_seasonal_cycle = gpcc_seasonal_cycle[:-1]
                                    array_storing_reanalysis[1] = gpcc_seasonal_cycle

                                pr_reference_array = np.mean(array_storing_reanalysis, axis = 0)

                                pr_ref_seasonal_cycle = [float('%.3f' % i) for i in pr_reference_array]

                    if variable == 'pr':

                        print "Ensemble"
                        ensemble_seasonal_cycle = np.asarray(ensemble_seasonal_cycle)
                        pr_ref_seasonal_cycle = np.asarray(pr_ref_seasonal_cycle)
                        print ensemble_seasonal_cycle
                        print pr_ref_seasonal_cycle
                        rmse_val = np.sqrt(((np.array(ensemble_seasonal_cycle[months]) - np.array(landfluxeval_seasonal_cycle[months])) ** 2).mean())
                        print rmse_val
                        pr_rmse = np.append(pr_rmse, rmse_val)
                        list_of_rmse_models_precip = np.append(list_of_rmse_models_precip, "Ensemble")

            if len(list_of_reanalysis) > 0:

                for i in np.arange(0, len(seasonal_cycle_reanalysis_array)):

                    reanalysis_number = i

                    reanalysis_id_1 = reanalysis_strings_for_plot[i]

                    reanalysis_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                    reanalysis_seasonal_cycle = [float('%.3f' % i) for i in reanalysis_seasonal_cycle]
                    reanalysis_seasonal_cycle = reanalysis_seasonal_cycle[:-1]

                    array_storing_reanalysis = np.zeros((2, 12))

                    for i in np.arange(0, len(seasonal_cycle_reanalysis_array)):

                        reanalysis_number = i

                        reanalysis_id_2 = reanalysis_strings_for_plot[i]

                        if variable == 'evaporation' and reanalysis_id_2 == 'LandFlux-EVAL':

                            landfluxeval_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                            landfluxeval_seasonal_cycle = [float('%.3f' % i) for i in landfluxeval_seasonal_cycle]
                            landfluxeval_seasonal_cycle = landfluxeval_seasonal_cycle[:-1]

                            if reanalysis_id_1 in ['LandFlux-EVAL', 'MODIS']:
                                pass
                            else:
                                if reanalysis_id_1 == 'GLEAM':
                                    reanalysis_id_1 = 'GLEAM/MSWEP'
                                print reanalysis_id_1

                                reanalysis_seasonal_cycle = np.asarray(reanalysis_seasonal_cycle)
                                landfluxeval_seasonal_cycle = np.asarray(landfluxeval_seasonal_cycle)

                                print reanalysis_seasonal_cycle[months]
                                print landfluxeval_seasonal_cycle[months]

                                rmse_val = np.sqrt(((np.array(reanalysis_seasonal_cycle[months]) - np.array(landfluxeval_seasonal_cycle[months])) ** 2).mean())
                                print rmse_val

                                evap_rmse = np.append(evap_rmse, rmse_val)
                                list_of_rmse_models_evap = np.append(list_of_rmse_models_evap, reanalysis_id_1)

                        if variable == 'pr':

                            if reanalysis_id_2 == 'CHIRPS':
                                chirps_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                                chirps_seasonal_cycle = [float('%.3f' % i) for i in chirps_seasonal_cycle]
                                chirps_seasonal_cycle = chirps_seasonal_cycle[:-1]
                                array_storing_reanalysis[0] = chirps_seasonal_cycle

                            if reanalysis_id_2 == 'GPCC':
                                gpcc_seasonal_cycle = seasonal_cycle_reanalysis_array[i]
                                gpcc_seasonal_cycle = [float('%.3f' % i) for i in gpcc_seasonal_cycle]
                                gpcc_seasonal_cycle = gpcc_seasonal_cycle[:-1]
                                array_storing_reanalysis[1] = gpcc_seasonal_cycle

                            pr_reference_array = np.mean(array_storing_reanalysis, axis = 0)

                            pr_ref_seasonal_cycle = [float('%.3f' % i) for i in pr_reference_array]

                    if variable == 'pr':

                        print reanalysis_id_1
                        reanalysis_seasonal_cycle = np.asarray(reanalysis_seasonal_cycle)
                        pr_ref_seasonal_cycle = np.asarray(pr_ref_seasonal_cycle)
                        print reanalysis_seasonal_cycle
                        print pr_ref_seasonal_cycle
                        rmse_val = np.sqrt(((np.array(reanalysis_seasonal_cycle[months]) - np.array(pr_ref_seasonal_cycle[months])) ** 2).mean())
                        print rmse_val

                        if reanalysis_id_1 in ['CHIRPS', 'GPCC']:
                            pass
                        else:
                            if reanalysis_id_1 == 'MSWEP':
                                reanalysis_id_1 = 'GLEAM/MSWEP'
                            list_of_rmse_models_precip = np.append(list_of_rmse_models_precip, reanalysis_id_1)
                            pr_rmse = np.append(pr_rmse, rmse_val)

            if variable == 'evaporation':
                list_of_rmse_models_evap = [str(i) for i in list_of_rmse_models_evap]
                print list_of_rmse_models_evap
                print evap_rmse

            if variable == 'pr':

                list_of_rmse_models_precip = [str(i) for i in list_of_rmse_models_precip]
                print list_of_rmse_models_precip
                print pr_rmse


        variable_number +=1

    if bar_plot == 'yes':

        print barchart_array_djf

        if 'DJF' in bar_times:
            print "hi1"
            plot_bar(barchart_array_djf, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'DJF', bar_width, variables_to_subtract)
        if 'MAM' in bar_times:
            print "hi2"
            plot_bar(barchart_array_mam, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'MAM', bar_width, variables_to_subtract)
        if 'JJA' in bar_times:
            print "hi3"
            print barchart_array_jja
            plot_bar(barchart_array_jja, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'JJA', bar_width, variables_to_subtract)
        if 'SON' in bar_times:
            print "hi4"
            plot_bar(barchart_array_son, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'SON', bar_width, variables_to_subtract)
        if 'Mar' in bar_times:
            print "hi4"
            plot_bar(barchart_array_mar, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'Mar', bar_width, variables_to_subtract)
        if 'Jul' in bar_times:
            print "hi4"
            plot_bar(barchart_array_jul, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'Jul', bar_width, variables_to_subtract)

    #print seasonal_cycle_reanalysis_multiple_variables_array

    if plot == 'yes':
        print "hi21"

        if number_of_variables > 1:

            print "hi5"
            model_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[0]
            reanalysis_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[1]
            print seasonal_cycle_models_multiple_variables_array
            print seasonal_cycle_ensemble_multiple_variables_array
            #print seasonal_cycle_ensemble_multiple_variables_std_array
            print seasonal_cycle_ensemble_multiple_variables_top_array
            print seasonal_cycle_ensemble_multiple_variables_bottom_array

            plot_seasonal_cycle(seasonal_cycle_models_multiple_variables_array, seasonal_cycle_ensemble_multiple_variables_array, seasonal_cycle_ensemble_multiple_variables_std_array, seasonal_cycle_ensemble_multiple_variables_top_array, seasonal_cycle_ensemble_multiple_variables_bottom_array, seasonal_cycle_dv_multiple_variables_array, seasonal_cycle_pv_multiple_variables_array, seasonal_cycle_reanalysis_multiple_variables_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim, y_tick_interval, number_of_variables, list_of_variables, variables_to_add, fill_between_lines, variables_to_subtract, legend_in_plot, number_of_models, list_of_models)

        if number_of_variables == 1 and subplot_multiple_variables == "yes":

            print "hi5"
            model_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[0]
            reanalysis_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[1]
            print seasonal_cycle_ensemble_multiple_variables_array
            print seasonal_cycle_ensemble_multiple_variables_top_array
            print seasonal_cycle_ensemble_multiple_variables_bottom_array

            if include_dynamic_prescribed == 'yes':
                print seasonal_cycle_dv_multiple_variables_array
                print seasonal_cycle_pv_multiple_variables_array

            plot_seasonal_cycle(seasonal_cycle_models_multiple_variables_array, seasonal_cycle_ensemble_multiple_variables_array, seasonal_cycle_ensemble_multiple_variables_std_array, seasonal_cycle_ensemble_multiple_variables_top_array, seasonal_cycle_ensemble_multiple_variables_bottom_array, seasonal_cycle_dv_multiple_variables_array, seasonal_cycle_pv_multiple_variables_array, seasonal_cycle_reanalysis_multiple_variables_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim, y_tick_interval, number_of_variables, list_of_variables, variables_to_add, fill_between_lines, variables_to_subtract, legend_in_plot, number_of_models, list_of_models)
