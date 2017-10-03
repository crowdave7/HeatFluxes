"""Import necessary modules for this code."""
import matplotlib
import matplotlib.pyplot as plt
import netCDF4
from netCDF4 import num2date
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
import os
matplotlib.use('Agg') 


def plot_bowen(list_of_models, model_type, list_of_months, season_name):
    """Import the data."""
    root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

    """Find the paths to the directories containing the model data"""
    directory_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in directories:
            path = os.path.join(root, i)
            for j in list_of_models:
                if j in path and model_type in path and ensemble in path:
                    directory_paths = np.append(directory_paths, path)

    """Find the model files for the sensible/latent heat fluxes and their absolute paths."""
    model_file_paths = []
    for i in directory_paths:
        files = os.listdir(i)
        for j in files:
            if "hfss" in j or "hfls" in j:
                model_file_path = os.path.join(i, j)
                model_file_paths = np.append(model_file_paths, model_file_path)
    model_file_paths.sort()
    #print model_file_paths

    """For the sensible heat flux files, extract the data."""

    for j in list_of_models:

        paths_for_this_model = [k for k in model_file_paths if j in k]
        sensible_path = [k for k in paths_for_this_model if 'hfss' in k]
        latent_path = [k for k in paths_for_this_model if 'hfls' in k]
        i = [sensible_path, latent_path]

        sensible_path = []
        latent_path = []

        if 'hfss' in i[0][0]:
            sensible_path = np.append(sensible_path, i[0][0])
        if 'hfls' in i[1][0]:
            latent_path = np.append(latent_path, i[1][0])

        sensible_data_model = extract_sensible_data(sensible_path[0], list_of_months)[0]

        latent_data_model = extract_latent_data(latent_path[0], list_of_months)
        print latent_data_model.shape

        evap_fraction = latent_data_model / (sensible_data_model + latent_data_model)

        """Begin to plot using Basemap."""

        model_id = extract_sensible_data(sensible_path[0], list_of_months)[1]
        fig = plt.figure()
        map = Basemap(llcrnrlon=-22, llcrnrlat=-22, urcrnrlon=62, urcrnrlat=12, projection='mill')

        """Draw coastlines, country lines, latitude and longitude lines"""
        map.drawcoastlines(linewidth=1.5)
        map.drawcountries(linewidth=1)
        map.drawparallels(np.arange(-50, 60, 10), labels=[1, 0, 0, 0], fontsize=10, linewidth=0.4)
        map.drawmeridians(np.arange(-40, 80, 20), labels=[0, 0, 0, 1], fontsize=10, linewidth=0.4)

        """Shift the grid to match up model data with basemap grid"""
        global_longitude = extract_sensible_data(sensible_path[0], list_of_months)[2]
        global_latitude = extract_sensible_data(sensible_path[0], list_of_months)[3]
        evap_fraction, global_longitude = shiftgrid(180., evap_fraction, global_longitude, start=False)
        grid = np.meshgrid(global_longitude, global_latitude)
        x, y = map(* grid)

        """Set contour levels, plot the data, and add a colour bar"""
        contour_levels = np.arange(0, 1.1, 0.1)
        contour_plot = map.contourf(x, y, evap_fraction, contour_levels, extend='both', cmap = 'coolwarm_r')
        colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        colour_bar.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        colour_bar.set_ticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)

        """Add title and label to the plot."""

        plt.title(model_id+' (AMIP) 1979-2008 ' ''+season_name+'', fontsize=10)
        colour_bar.set_label(("Evaporative Fraction"), fontsize=10)

        """Save the figure, and close the plot"""
        fig.savefig("Evap_Frac_"+season_name+"_"+model_id+".png")
        plt.close()
        print model_id+" plot done"

def extract_sensible_data(file_path, list_of_months):
    data = netCDF4.Dataset(file_path)
    model_id = data.model_id

    """Print the times."""
    time_data = num2date(data.variables['time'][:], units = data.variables['time'].units, calendar = data.variables['time'].calendar)

    #print len(time_data)

    """Find the relevant indices to slice the time data."""
    indices_list_time = []
    for index, i in enumerate(time_data):
            if i.month in (list_of_months):
                indices_list_time = np.append(indices_list_time, (index))
    indices_list_time = [int(i) for i in indices_list_time]
    #print indices_list_time

    if 1 and 12 in list_of_months:
        indices_list_time = [(i+12) if (1 <= time_data[i].month <= 6) else i for i in indices_list_time]
        indices_list_time = [i for i in indices_list_time if i < (len(time_data)-1)]

    print indices_list_time

    """Compute the seasonal mean sensible heat flux."""
    sensible_data_model = data.variables["hfss"][indices_list_time].mean(axis = 0)

    global_longitude = data.variables['lon'][:]
    global_latitude = data.variables['lat'][:]

    return sensible_data_model, model_id, global_longitude, global_latitude

def extract_latent_data(file_path, list_of_months):
    data = netCDF4.Dataset(file_path)

    """Print the times."""
    time_data = num2date(data.variables['time'][:], units = data.variables['time'].units, calendar = data.variables['time'].calendar)

    """Find the relevant indices to slice the time data."""
    indices_list_time = []
    for index, i in enumerate(time_data):
            if i.month in (list_of_months):
                indices_list_time = np.append(indices_list_time, (index))
    indices_list_time = [int(i) for i in indices_list_time]

    if 1 and 12 in list_of_months:
        indices_list_time = [(i+12) if (1 <= time_data[i].month <= 6) else i for i in indices_list_time]
        indices_list_time = [i for i in indices_list_time if i < (len(time_data)-1)]

    """Compute the seasonal mean latent heat flux."""
    latent_data_model = data.variables["hfls"][indices_list_time].mean(axis = 0)

    return latent_data_model

plot_bowen(["GFDL-CM3/"], "amip", [1,2,12], "DJF")
