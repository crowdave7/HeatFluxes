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


def plot_map(list_of_models, model_type, variable, list_of_months, season_name):
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

    """Find the model files and their absolute paths."""
    model_file_paths = []
    for i in directory_paths:
        files = os.listdir(i)
        for j in files:
            if variable in j:
                model_file_path = os.path.join(i, j)
                model_file_paths = np.append(model_file_paths, model_file_path)
    model_file_paths.sort()

    """For each model file, plot the map."""
    for i in model_file_paths:

        print i

        """Extract the data, model ID, variable name and variable units."""
        data = netCDF4.Dataset(i)
        model_id = data.model_id
        variable_name = data.variables[variable].long_name
        variable_units = data.variables[variable].units

        """Print a list of the variables."""
        # print data.variables

        """Find the dimensions and their sizes."""
        # for i in data.dimensions:
            # print data.dimensions[str(i)].name
            # print data.dimensions[str(i)].size

        """Dimensions sensible heat flux lon (144), bnds (2), lat (90), time (360)"""

        """Find the variables and their sizes."""
        # for i in data.variables:
            # print data.variables[str(i)].name
            # print data.variables[str(i)].size

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

        """Compute the seasonal mean sensible heat flux."""
        mean_field = data.variables[variable][indices_list_time].mean(axis = 0)

        """Begin to plot using Basemap."""
        fig = plt.figure()
        #map = Basemap(llcrnrlon=5, llcrnrlat=-10, urcrnrlon=35, urcrnrlat=5, projection='mill')
        map = Basemap(llcrnrlon=-22, llcrnrlat=-22, urcrnrlon=62, urcrnrlat=12, projection='mill')
        """Draw coastlines, country lines, latitude and longitude lines"""
        map.drawcoastlines(linewidth=1.5)
        map.drawcountries(linewidth=1)
        map.drawparallels(np.arange(-50, 60, 10), labels=[1, 0, 0, 0], fontsize=10, linewidth=0.4)
        map.drawmeridians(np.arange(-40, 80, 20), labels=[0, 0, 0, 1], fontsize=10, linewidth=0.4)

        """Shift the grid to match up model data with basemap grid"""
        longitude = data.variables['lon'][:]
        latitude = data.variables['lat'][:]
        mean_field, longitude = shiftgrid(180., mean_field, longitude, start=False)
        grid = np.meshgrid(longitude, latitude)
        x, y = map(* grid)

        """Set contour levels, plot the data, and add a colour bar"""
        contour_levels = np.arange(0, 165, 15)
        contour_plot = map.contourf(x, y, mean_field, contour_levels, extend='both', cmap = 'coolwarm')
        colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        colour_bar.set_ticks([0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150])
        colour_bar.set_ticklabels([0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150])
        colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)

        #SH 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 0, 110, 10
        #LH 0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150 0, 165, 15

        """Add title and label to the plot."""

        plt.title(model_id+' (AMIP) 1979-2008 ' ''+season_name+'', fontsize=10)
        colour_bar.set_label(variable_name+" ("+variable_units+")", fontsize=10)

        """Save the figure, and close the plot"""
        fig.savefig(variable+"_"+season_name+"_"+model_id+".png")
        plt.close()
        print model_id+" plot done"

plot_map(["GFDL-CM3/"], "amip", "hfls", [1,2,12], "DJF")
