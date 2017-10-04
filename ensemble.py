import matplotlib
import matplotlib.pyplot as plt
import netCDF4
from netCDF4 import num2date
from netCDF4 import date2num
import numpy as np
from mpl_toolkits.basemap import Basemap, maskoceans
from mpl_toolkits.basemap import shiftgrid
import os
import iris
import iris.analysis
import iris.analysis.cartography
import copy
matplotlib.use('Agg')


def plot_ensemble_mean_map(variable, list_of_months, list_of_models):

    root_directory = "/ouce-home/students/kebl4396/Paper1/regridded"

    list_of_cubes = []

    model_file_paths = []

    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            #print path
            for j in list_of_models:
                #Slice slash off inputted models
                j = j[:-1]
                if j in path and variable in path:
                    #print path
                    model_file_paths = np.append(model_file_paths, path)

    #print model_file_paths

    for i in model_file_paths:

        if 'hfls' in i:
            name = 'surface_upward_latent_heat_flux'
            #print name
            #print i
            data = iris.load_cube(i, name)
            times = extract_times(netCDF4.Dataset(i), list_of_months)

        if 'hfss' in i:
            name = 'surface_upward_sensible_heat_flux'
            data = iris.load_cube(i, name)
            times = extract_times(netCDF4.Dataset(i), list_of_months)

        data = data.extract(iris.Constraint(time = times))

        data_unmasked = data.collapsed('time', iris.analysis.MEAN)

        #longitude = data_unmasked.coord('longitude').points

        #latitude = data_unmasked.coord('latitude').points

        #print longitude
        #print latitude

        #print data_unmasked

        data = data_unmasked

        #print data

        list_of_cubes = np.append(list_of_cubes, data)

    ensemble_mean = sum(list_of_cubes)/len(list_of_cubes)
    print ensemble_mean

def extract_times(data, list_of_months):

    #Print the times.

    time_data = num2date(data.variables['time'][:], units = data.variables['time'].units)

    list_time = np.array(time_data).tolist()

    #print list_time

    indices_list_time = []

    for index, i in enumerate(list_time):

        for elem in list_of_months:

           if i.month == elem:

                indices_list_time = np.append(indices_list_time, (index))

    indices_list_time = [int(i) for i in indices_list_time]

    time_data = time_data[indices_list_time]

    time_data = date2num(time_data, data.variables['time'].units)

    np.set_printoptions(suppress = True)

    return time_data


plot_ensemble_mean_map("hfls", [1,2,3,4,5,6,7,8,9,10,11,12], ["GFDL-CM3/", "HadGEM2-A/"])
