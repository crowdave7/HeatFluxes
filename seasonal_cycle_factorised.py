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

def model_file_paths_func(list_of_models, model_type, variable):

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

    if variable == 'nrad':
        print "hi9"
        """Find the paths to the files containing the model data"""
        model_file_paths = []
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

    if variable not in ["nrad", "evap_fraction"]:

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

        model_file_paths = model_file_paths_sorted

    print "hi2"
    print model_file_paths
    return model_file_paths

def slicing(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        if variable != "treeFrac":
            cubelist[i] = cubelist[i].extract(time_range)
            time_points = cubelist[i].coord('time').points
            times = cubelist[i].coord('time').units.num2date(time_points)
        if variable == 'nrad':
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
    if variable == 'nrad':
        model_id = cubelist[i].long_name
    else:
        model_id = cubelist[i].attributes['model_id']

    if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg' or variable == 'mrros' or variable == 'mrro' or variable == 'prveg':
        cubelist[i] = iris.analysis.maths.multiply(cubelist[i], 86400)

    """if variable is evaporation, divide by 28"""
    if variable == 'evaporation':
        if unit_plot == 'mm day-1':
            cubelist[i] = iris.analysis.maths.divide(cubelist[i], 28)

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
            cubelist_ensemble[i] = cubelist_ensemble[i].extract(time_range)
            time_points = cubelist_ensemble[i].coord('time').points
            times = cubelist_ensemble[i].coord('time').units.num2date(time_points)
        if variable == 'nrad':
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
    if variable == 'nrad':
        model_id = cubelist_ensemble[i].long_name
    else:
        model_id = cubelist_ensemble[i].attributes['model_id']

    if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg' or variable == 'mrros' or variable == 'mrro' or variable == 'prveg':
        cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 86400)

    """if variable is evaporation, divide by 28"""
    if variable == 'evaporation':
        if unit_plot == "mm day-1":
            cubelist_ensemble[i] = iris.analysis.maths.divide(cubelist_ensemble[i], 28)

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
                cubelist_reanalysis[i] = cubelist_reanalysis[i].extract(time_range)
            if reanalysis_id == "era5":
                cubelist_reanalysis[i] = cubelist_reanalysis[i].extract(iris.Constraint(time=lambda cell: 2010 <= cell.point.year <= 2017))
            time_points = cubelist_reanalysis[i].coord('time').points
            times = cubelist_reanalysis[i].coord('time').units.num2date(time_points)
        if variable != "treeFrac":
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
        cubelist_reanalysis[i].long_name = "GLEAM-LE"
        cubelist_reanalysis[i].rename("GLEAM-LE")
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


def plot_seasonal_cycle(seasonal_cycle_models_array, seasonal_cycle_ensemble_array, seasonal_cycle_reanalysis_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim, y_tick_interval, lower_ylim_right, upper_ylim_right, y_tick_interval_right, number_of_variables, list_of_variables, variables_to_add, fill_between_lines, variables_to_subtract, legend_in_plot):

    """Before looping through all the models, set up the figure to plot to."""

    print "Setting up figure"
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    objects = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan')
    x_pos = np.arange(len(objects))
    plt.xticks(x_pos, objects)
    ax1.set_xlim([0, 11])
    ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
    ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')
    ax1.patch.set_visible(False)

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

                seasonal_cycle_model = seasonal_cycle_models_array[i]
                line_colour = model_line_colours[i]
                model_id = model_strings_for_plot[i]
                print model_id

                ax1.plot(x_pos, seasonal_cycle_model, zorder=1, linestyle='-', color=line_colour, label = str(model_id))
                handles, labels = ax1.get_legend_handles_labels()
                handles[-1].set_linestyle("-")
                legend = plt.legend(handles, labels, loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9, handlelength=2.5)

            """Plot model ensemble mean seasonal cycle."""

        if len(ensemble_string_for_plot) == 0:
            pass

        else:

            ensemble_mean_seasonal_cycle_models = seasonal_cycle_ensemble_array
            ensemble_string = ensemble_string_for_plot[0]
            print ensemble_string

            ax1.plot(x_pos, ensemble_mean_seasonal_cycle_models, zorder=2, linestyle='-', linewidth=4.0, color='black', label = str(ensemble_string))
            handles, labels = ax1.get_legend_handles_labels()
            handles[-1].set_linestyle("-")
            legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

            """Plot reanalysis seasonal cycles."""

            if len(reanalysis_strings_for_plot) == 0:
                pass
            else:

                for i in np.arange(0, len(seasonal_cycle_reanalysis_array)):

                    seasonal_cycle_reanalysis = seasonal_cycle_reanalysis_array[i]
                    line_colour = reanalysis_line_colours[i]
                    reanalysis_id = reanalysis_strings_for_plot[i]
                    print reanalysis_id

                    ax1.plot(x_pos, seasonal_cycle_reanalysis, zorder=1, linestyle='--', color = line_colour, label = str(reanalysis_id))
                    handles, labels = ax1.get_legend_handles_labels()
                    handles[-1].set_linestyle("--")
                    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

        print lower_y_lim
        print upper_y_lim+y_tick_interval
        print y_tick_interval
        print np.arange(lower_y_lim, upper_y_lim+y_tick_interval, y_tick_interval)
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
            plt.ylabel('Total Runoff Flux (mm $\mathregular{day^{-1}}$)')


        """Save the figure."""
        print "Saving figure"
        fig.savefig("Seasonal_Cycle_"+variable+".png", bbox_extra_artists=(legend,), bbox_inches='tight')
        print "Plot done."

    if number_of_variables > 1:
        print "hi8"
        print "Figure set up"
        print len(variables_to_add)
        # for i in variables_to_add:
        #     for j in np.arange(0, len(i)):
        #         print i[j]

        """Load the data from the model file paths into a cube. Constrain the input years"""
        """Print the model ID, length of time dimension, and first and last model dates."""

        """Plot model seasonal cycles."""

        if len(reanalysis_strings_for_plot) > 0:

            variable_number = 0
            line_number = 0
            line_count = 0

            fill_between_lines_data = seasonal_cycle_to_add = np.zeros((2, 12))

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
                    seasonal_cycle_to_add = np.zeros((len(i), 12))
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
                    seasonal_cycle_to_subtract = np.zeros((len(i), 12))
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

        if len(ensemble_string_for_plot) > 0:

            variable_number = 0
            line_number = 0
            line_count = 0

            fill_between_lines_data = seasonal_cycle_to_add = np.zeros((2, 12))

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
                    color = 'red'
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
                    seasonal_cycle_to_add = np.zeros((len(i), 12))
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
                    seasonal_cycle_to_subtract = np.zeros((len(i), 12))
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
            plt.ylabel('Precipitation (mm $\mathregular{day^{-1}}$)')

        """Save the figure."""
        print "Saving figure"
        fig.savefig("Seasonal_Cycle_"+variable+".png", bbox_inches='tight')
        print "Plot done."


def plot_bar(barchart_array, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, bar_season, bar_width, variables_to_subtract):

    if number_of_variables == 3:

        fig = plt.figure()

        variable_data_1 = barchart_array[0,:]
        variable_data_2 = barchart_array[1,:]
        variable_data_3 = barchart_array[2,:]

        if ensemble == 'yes':
            number_of_models = number_of_models+1

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

        if ensemble == 'yes':
            model_strings_for_plot = np.append(model_strings_for_plot, "Ensemble")
        if number_of_reanalysis > 0:
            model_strings_for_plot = np.append(model_strings_for_plot, reanalysis_strings_for_plot)

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

            list_of_variables.insert(3, string)

            print list_of_variables
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
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = []

    list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "BNU-ESM/"]
    #list_of_models = ["bcc-csm1-1/"]
    #list_of_models = []
    list_of_reanalysis = []
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #list_of_variables = ["evaporation", "tran", "evspsblveg", "evspsblsoi"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['pr']
    #list_of_variables = ['pr']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    list_of_variables = ['pr', 'evaporation', 'mrro']

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    variables_to_subtract = [['pr', 'mrro']]

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    fill_between_lines = [1, 4]
    #variables_to_add = []

    lower_lat = -10
    upper_lat = 5
    lower_lon = 5
    upper_lon = 35
    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 8.0
    y_tick_interval = 1.0
    lower_ylim_right = 0
    upper_ylim_right = 100
    y_tick_interval_right = 10

    cmap = 'rainbow'
    ensemble = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_seasons = ['DJF']
    bar_y_axis_title = 'Precipitation (mm $\mathregular{day^{-1}}$)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'


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

    if ensemble == 'yes':

        barchart_array_djf = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_mam = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_jja = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))
        barchart_array_son = np.zeros((len(list_of_variables), (len(list_of_models)+len(list_of_reanalysis)+1)))

    """For each variable, calculate the model, ensemble and reanalysis seasonal cycle arrays. """

    """Define first variable number."""

    variable_number = 0

    seasonal_cycle_models_multiple_variables_array = []
    seasonal_cycle_ensemble_multiple_variables_array = np.zeros((len(list_of_variables), 12))
    seasonal_cycle_reanalysis_multiple_variables_array = np.zeros((len(list_of_variables), 12))

    for i in list_of_variables:

        seasonal_cycle_models_array = []
        seasonal_cycle_ensemble_array = []
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

                """Extract the regridded model file paths for the ensemble mean."""
                model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

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
                    p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_ensemble = array.values()

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


                ensemble_string_for_plot = np.append(ensemble_string_for_plot, "Ensemble")

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

        if plot == 'yes':

            if number_of_variables == 1:

                model_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[0]
                reanalysis_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[1]
                plot_seasonal_cycle(seasonal_cycle_models_array, seasonal_cycle_ensemble_array, seasonal_cycle_reanalysis_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim, y_tick_interval, lower_ylim_right, upper_ylim_right, y_tick_interval_right, number_of_variables, list_of_variables, variables_to_add, fill_between_lines, variables_to_subtract, legend_in_plot)
                print('')

    #------------------------------------------------------------------------------------------------------------------------------------------

        """PRINT VALUES FOR SEASONAL CYCLE."""

        djf_values = []
        mam_values = []
        jja_values = []
        son_values = []
        list_all_datasets = []

        if len(list_of_models) > 0:
            list_all_datasets = np.append(list_all_datasets, list_of_models)
        if ensemble == "yes":
            list_all_datasets = np.append(list_all_datasets, "Ensemble")
        if len(list_of_reanalysis) > 0:
            list_all_datasets = np.append(list_all_datasets, list_of_reanalysis)

        print "hi3"
        print list_of_reanalysis
        print list_all_datasets

        if len(list_of_models) == 0:
            model_number = -1

        if len(list_of_models) > 0:

            for i in np.arange(0, len(seasonal_cycle_models_array)):

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

                print "%.3f, %.3f, %.3f, %.3f" % (djf_value, mam_value, jja_value, son_value)

                print('')

                djf_values = np.append(djf_values, djf_value)
                mam_values = np.append(mam_values, mam_value)
                jja_values = np.append(jja_values, jja_value)
                son_values = np.append(son_values, son_value)

                barchart_array_djf[variable_number, model_number] = djf_value
                barchart_array_mam[variable_number, model_number] = mam_value
                barchart_array_jja[variable_number, model_number] = jja_value
                barchart_array_son[variable_number, model_number] = son_value

                print variable_number
                print model_number

        if ensemble == 'yes':

            print "Ensemble"
            ensemble_seasonal_cycle = [float('%.3f' % i) for i in seasonal_cycle_ensemble_array]
            print ensemble_seasonal_cycle

            djf_value = (ensemble_seasonal_cycle[11] + ensemble_seasonal_cycle[0] + ensemble_seasonal_cycle[1])/float(3.0)
            mam_value = (ensemble_seasonal_cycle[2] + ensemble_seasonal_cycle[3] + ensemble_seasonal_cycle[4])/float(3.0)
            jja_value = (ensemble_seasonal_cycle[5] + ensemble_seasonal_cycle[6] + ensemble_seasonal_cycle[7])/float(3.0)
            son_value = (ensemble_seasonal_cycle[8] + ensemble_seasonal_cycle[9] + ensemble_seasonal_cycle[10])/float(3.0)

            print "%.3f, %.3f, %.3f, %.3f" % (djf_value, mam_value, jja_value, son_value)

            #print('')

            djf_values = np.append(djf_values, djf_value)
            mam_values = np.append(mam_values, mam_value)
            jja_values = np.append(jja_values, jja_value)
            son_values = np.append(son_values, son_value)

            model_number = model_number+1

            barchart_array_djf[variable_number, model_number] = djf_value
            barchart_array_mam[variable_number, model_number] = mam_value
            barchart_array_jja[variable_number, model_number] = jja_value
            barchart_array_son[variable_number, model_number] = son_value

            print variable_number

            seasonal_cycle_ensemble_multiple_variables_array[variable_number] = ensemble_seasonal_cycle

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

                print "%.3f, %.3f, %.3f, %.3f" % (djf_value, mam_value, jja_value, son_value)

                print('')

                djf_values = np.append(djf_values, djf_value)
                mam_values = np.append(mam_values, mam_value)
                jja_values = np.append(jja_values, jja_value)
                son_values = np.append(son_values, son_value)

                reanalysis_number = reanalysis_number + 1 + model_number

                barchart_array_djf[variable_number, reanalysis_number] = djf_value
                barchart_array_mam[variable_number, reanalysis_number] = mam_value
                barchart_array_jja[variable_number, reanalysis_number] = jja_value
                barchart_array_son[variable_number, reanalysis_number] = son_value

                seasonal_cycle_reanalysis_multiple_variables_array[variable_number] = reanalysis_seasonal_cycle

        list_all_datasets_input = [str(i) for i in list_all_datasets]

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


        """
        print mam_values
        print jja_values
        print son_values
        print list_all_datasets
        """

        variable_number +=1

    if bar_plot == 'yes':

        print barchart_array_djf

        if 'DJF' in bar_seasons:
            print "hi1"
            plot_bar(barchart_array_djf, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'DJF', bar_width, variables_to_subtract)
        if 'MAM' in bar_seasons:
            print "hi2"
            plot_bar(barchart_array_mam, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'MAM', bar_width, variables_to_subtract)
        if 'JJA' in bar_seasons:
            print "hi3"
            print barchart_array_jja
            plot_bar(barchart_array_jja, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'JJA', bar_width, variables_to_subtract)
        if 'SON' in bar_seasons:
            print "hi4"
            plot_bar(barchart_array_son, number_of_models, number_of_reanalysis, number_of_variables, list_of_variables, model_strings_for_plot, reanalysis_strings_for_plot, ensemble, bar_y_axis_title, bar_colours, lower_y_lim, upper_y_lim, 'SON', bar_width, variables_to_subtract)

    #print seasonal_cycle_ensemble_multiple_variables_array
    print seasonal_cycle_reanalysis_multiple_variables_array
    if plot == 'yes':

        if number_of_variables > 1:

            print "hi5"
            model_line_colours = []
            reanalysis_line_colours = []
            plot_seasonal_cycle(seasonal_cycle_models_multiple_variables_array, seasonal_cycle_ensemble_multiple_variables_array, seasonal_cycle_reanalysis_multiple_variables_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim, y_tick_interval, lower_ylim_right, upper_ylim_right, y_tick_interval_right, number_of_variables, list_of_variables, variables_to_add, fill_between_lines, variables_to_subtract, legend_in_plot)
