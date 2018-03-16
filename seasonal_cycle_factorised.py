#!/usr/bin/env python

"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans, interp
import multiprocessing
import numpy as np
import os
import time
matplotlib.use('Agg')


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

    """If variable is mrso, distinguish between mrso and mrsos to find model files."""
    if variable == 'mrso':
        variable = 'mrso_'

    """If variable is mrro, distinguish between mrro and mrros to find model files."""
    if variable == 'mrro':
        variable = 'mrro_'

    if variable == 'evapotranspiration':
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
        variable = 'evapotranspiration'

    return model_file_paths


def find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory):

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

        model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    if variable == "evap_fraction":

        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                for j in list_of_models:
                    if j in path and model_type in path and ensemble in path and ('hfss' in path or 'hfls' in path):
                        model_file_paths = np.append(model_file_paths, path)

        model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    if variable not in ["nrad", "evap_fraction"]:

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

def run_cube_year_slicing_change_units(i, array):

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

        print i
        print data

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

        print i
        print data

        """Append the data to the array outside the loop to produce the data for the seasonal cycle."""
        array[i] = data

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

        """If variable is mrso, distinguish between mrso and mrsos to find model files."""
        if variable == 'mrso':
            variable = 'mrso_'

        """If variable is mrro, distinguish between mrro and mrros to find model files."""
        if variable == 'mrro':
            variable = 'mrro_'

        if variable == 'evapotranspiration':
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

def run_cube_year_slicing_change_units_ensemble(i, array):

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

    """If variable is mrso, distinguish between mrso and mrsos to find model files."""
    if variable == 'mrso':
        variable = 'mrso_'

    """If variable is mrro, distinguish between mrro and mrros to find model files."""
    if variable == 'mrro':
        variable = 'mrro_'

    if variable == 'evapotranspiration':
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


def run_cube_year_slicing_change_units_reanalysis(i, array):

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

    """Slice the regridded cube down to the African domain."""
    cubelist_reanalysis[i] = cubelist_reanalysis[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

    """if variable is evapotranspiration, divide by 28"""
    if variable == 'evapotranspiration':
        cubelist_reanalysis[i] = iris.analysis.maths.divide(cubelist_reanalysis[i], 28)

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


def plot_seasonal_cycle(seasonal_cycle_models_array, seasonal_cycle_ensemble_array, seasonal_cycle_reanalysis_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim):

    """Before looping through all the models, set up the figure to plot to."""

    print "Setting up figure"
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    objects = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec')
    x_pos = np.arange(len(objects))
    plt.xticks(x_pos, objects)
    ax1.set_xlim([0, 11])
    ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
    ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')
    ax1.patch.set_visible(False)

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

        ensemble_mean_seasonal_cycle_models = seasonal_cycle_ensemble_array
        ensemble_string = ensemble_string_for_plot
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

    """Add lower and upper y limits."""
    plt.ylim(lower_y_lim, upper_y_lim)

    if variable == 'pr':
        plt.ylabel('Precipitation (mm $\mathregular{day^{-1}}$)')
    if variable == 'hfls':
        plt.ylabel('Surface Upward Latent Heat Flux (W $\mathregular{m^{-2}}$)')
    if variable == 'evapotranspiration':
        plt.ylabel('Evapotranspiration (mm $\mathregular{day^{-1}}$)')
    if variable == 'hfss':
        plt.ylabel('Surface Upward Sensible Heat Flux (W $\mathregular{m^{-2}}$)')
    if variable == 'evap_fraction':
        plt.ylabel('Evaporative Fraction')
    if variable == 'nrad':
        plt.ylabel('Surface Net Radiation (W $\mathregular{m^{-2}}$)')
    if variable == 'mrsos':
        plt.ylabel('Soil Moisture Content of Upper Layer (mm)')
    if variable == 'mrso':
        plt.ylabel('Soil Moisture Content (mm)')
    if variable == 'tran':
        plt.ylabel('Transpiration (mm $\mathregular{day^{-1}}$)')
    if variable == 'evspsbl':
        plt.ylabel('Evaporation (mm $\mathregular{day^{-1}}$)')
    if variable == 'evspsblsoi':
        plt.ylabel('Bare Soil Evaporation (mm $\mathregular{day^{-1}}$)')
    if variable == 'evspsblveg':
        plt.ylabel('Evaporation from Canopy (mm $\mathregular{day^{-1}}$)')
    if variable == 'prveg':
        plt.ylabel('Precipitation Intercepted by Canopy (mm $\mathregular{day^{-1}}$)')
    if variable == 'mrros':
        plt.ylabel('Surface Runoff Flux (mm $\mathregular{day^{-1}}$)')
    if variable == 'lai':
        plt.ylabel('Leaf Area Index')
    if variable == 'mrro':
        plt.ylabel('Runoff Flux (mm $\mathregular{day^{-1}}$)')

    """Save the figure."""
    print "Saving figure"
    fig.savefig("Seasonal_Cycle_"+variable+".png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)
    print "Plot done."

if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "GISS-E2-R/"]
    model_type = "amip"
    list_of_reanalysis = ["cfsr", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = []
    variable = "evapotranspiration"
    lower_lat = -10
    upper_lat = 5
    lower_lon = 5
    upper_lon = 35
    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 1
    upper_y_lim = 5
    cmap = 'rainbow'
    ensemble = 'yes'
    plot = 'no'

    # ------------------------------------------------------------------------------------------------------------------------------------------

    time_range = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)

    seasonal_cycle_models_array = []
    seasonal_cycle_ensemble_array = []
    seasonal_cycle_reanalysis_array = []
    model_strings_for_plot = []
    ensemble_string_for_plot = "Ensemble"
    reanalysis_strings_for_plot = []

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # MODELS AND ENSEMBLE MEAN

    if len(list_of_models) == 0:
        pass
    else:

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

        """Slice each cube by the time range."""
        manager = multiprocessing.Manager()
        array = manager.dict()
        jobs = []
        for i in np.arange(0, len(cubelist)):
            p = multiprocessing.Process(target=run_cube_year_slicing_change_units, args=(i, array))
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
                seasonal_cycle_models_array.append(model_seasonal_cycle)

        if ensemble == 'yes':

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
                p = multiprocessing.Process(target=run_cube_year_slicing_change_units_ensemble, args=(i, array))
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

                model_id = cubelist[i].long_name

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
                    seasonal_cycle_array_for_ensemble.append(cube_seasonal_cycle_for_ensemble)

            seasonal_cycle_ensemble_array = np.mean(seasonal_cycle_array_for_ensemble, axis = 0)

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
            p = multiprocessing.Process(target=run_cube_year_slicing_change_units_reanalysis, args=(i, array))
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
                reanalysis_strings_for_plot = np.append(reanalysis_strings_for_plot, reanalysis_id)
                seasonal_cycle_reanalysis_array.append(cube_seasonal_cycle_reanalysis)

    if len(list_of_models) > 0:
        print np.shape(seasonal_cycle_models_array)
    if ensemble == 'yes':
        print np.shape(seasonal_cycle_ensemble_array)
    if len(list_of_reanalysis) > 0:
        print np.shape(seasonal_cycle_reanalysis_array)
    if len(list_of_models) > 0:
        print model_strings_for_plot
    if len(list_of_reanalysis) > 0:
        print reanalysis_strings_for_plot

    number_of_models = len(list_of_models)
    number_of_reanalysis = len(list_of_reanalysis)

    if plot == 'yes':
        model_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[0]
        reanalysis_line_colours = line_colours(number_of_models, number_of_reanalysis, cmap)[1]
        plot_seasonal_cycle(seasonal_cycle_models_array, seasonal_cycle_ensemble_array, seasonal_cycle_reanalysis_array, model_strings_for_plot, ensemble_string_for_plot, reanalysis_strings_for_plot, model_line_colours, reanalysis_line_colours, lower_y_lim, upper_y_lim)

    if len(list_of_models) > 0:

        for i in np.arange(0, len(seasonal_cycle_models_array)):

            model_id = model_strings_for_plot[i]
            print model_id
            print seasonal_cycle_models_array[i]

        if ensemble == 'yes':
            print "Ensemble"
            print seasonal_cycle_ensemble_array

    if len(list_of_reanalysis) > 0:

        for i in np.arange(0, len(seasonal_cycle_reanalysis_array)):

            reanalysis_id = reanalysis_strings_for_plot[i]
            print reanalysis_id
            print seasonal_cycle_reanalysis_array[i]
