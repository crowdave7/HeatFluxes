"""Import necessary modules for this code."""
import copy
import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap, maskoceans
import netCDF4
from netCDF4 import num2date
from netCDF4 import date2num
import numpy as np
import os
import seasonal_cycle_ensemble
matplotlib.use('Agg')


def seasonal_cycle(list_of_models, model_type, variable, lower_lat, upper_lat, lower_lon, upper_lon):
    """Take the input variables, and find the paths to the relevant regridded model files."""
    """Plot these up as a seasonal cycle. Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

    if variable != "evap_fraction":

        """Find the paths to the files containing the model data"""
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

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    """Before looping through all the models, set up the figure to plot to."""

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    objects = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec')
    x_pos = np.arange(len(objects))
    plt.xticks(x_pos, objects)
    ax1.set_xlim([0, 11])
    ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
    ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')
    ax1.patch.set_visible(False)
    plt.title('Central Africa, AMIP, 1979-2008', fontsize = 10)

    """Load the data from the model file paths into a cube. Constrain the input years"""
    """Print the model ID, length of time dimension, and first and last model dates."""

    """Define the list of colours for each model."""
    cmap = plt.get_cmap('rainbow')
    colours = [cmap(i) for i in np.linspace(0, 1, len(list_of_models))]

    """Define a model number to begin with."""
    model_number = 0

    """For each model,"""
    for j in list_of_models:

        """Select the line colour and add one to the line count variable."""
        line_colour = colours[model_number]

        """Extract the sensible and latent paths."""

        paths_for_this_model = [k for k in model_file_paths if j in k]
        #print paths_for_this_model
        latent_path = [k for k in paths_for_this_model if 'hfls' in k]
        sensible_path = [k for k in paths_for_this_model if 'hfss' in k]

        """Extract the model ID."""
        data = netCDF4.Dataset(paths_for_this_model[0])
        model_id = data.model_id

        """Load the cube for each variable, constrain the years and extract the seasonal cycle array."""
        """Append the seasonal cycle array to the ensemble mean array outside the loop."""

        if variable == 'hfls':
            data_cube = iris.load_cube(latent_path, 'surface_upward_latent_heat_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, latent_path[0], lower_lat, upper_lat, lower_lon, upper_lon)

        if variable == 'hfss':
            data_cube = iris.load_cube(sensible_path, 'surface_upward_sensible_heat_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, sensible_path[0], lower_lat, upper_lat, lower_lon, upper_lon)

        if variable == 'evap_fraction':
            cube_latent = iris.load_cube(latent_path, 'surface_upward_latent_heat_flux')
            cube_latent = constrain_year(cube_latent, time_range)
            cube_sensible = iris.load_cube(sensible_path, 'surface_upward_sensible_heat_flux')
            cube_sensible = constrain_year(cube_sensible, time_range)

            latent_seasonal_cycle_array = extract_seasonal_cycle_data(cube_latent, latent_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            sensible_seasonal_cycle_array = extract_seasonal_cycle_data(cube_sensible, sensible_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            data_array = latent_seasonal_cycle_array / (latent_seasonal_cycle_array + sensible_seasonal_cycle_array)

        """Add the seasonal cycle to the plot."""

        ax1.plot(x_pos, data_array, zorder=1, linestyle='-', color=line_colour, label = str(model_id))
        handles, labels = ax1.get_legend_handles_labels()
        handles = [copy.copy(ha) for ha in handles]
        [ha.set_linestyle("-") for ha in handles]
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1, 1))
        if variable == 'hfls':
            plt.ylabel('Latent Heat Flux (W m-2)')
            plt.ylim(0, 160)
        if variable == 'hfss':
            plt.ylabel('Sensible Heat Flux (W m-2)')
            plt.ylim(0, 160)
        if variable == 'evap_fraction':
            plt.ylabel('Evaporative Fraction')
            plt.ylim(0, 1)

        """Add 1 to the model number to loop through the next model."""
        model_number +=1

    """Take the mean seasonal cycle across the models."""

    ensemble_mean_array = seasonal_cycle_ensemble.seasonal_cycle_ensemble(list_of_models, model_type, variable, lower_lat, upper_lat, lower_lon, upper_lon)
    print ensemble_mean_array

    """Add the ensemble mean seasonal cycle to the plot."""

    ax1.plot(x_pos, ensemble_mean_array, zorder=1, linestyle=':', color='black', label = "Ensemble")
    handles, labels = ax1.get_legend_handles_labels()
    handles = [copy.copy(ha) for ha in handles]
    [ha.set_linestyle("-") for ha in handles]
    legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1, 1))

    """Save the figure."""
    fig.savefig("Seasonal_Cycle_"+variable+"_test.png", bbox_extra_artists=(legend,), bbox_inches='tight')


def extract_seasonal_cycle_data(input_cube, path, lower_lat, upper_lat, lower_lon, upper_lon):
    """Extract seasonal cycle data from an iris cube given the cube, its path, and the chosen lats/lons."""
    data_array = []

    """For each month,"""
    for i in np.arange(1, 13, 1):

        """Load the dataset using the input path, and load the times to constrain the data by month."""
        times = constrain_month(netCDF4.Dataset(path), [i])

        """Load in the cube."""
        data = input_cube

        """Constrain the times in the cube."""
        data = data.extract(iris.Constraint(time = times))

        """Constrain the latitudes and longitudes of the data."""
        data = data.intersection(latitude=(lower_lat, upper_lat), longitude=(lower_lon, upper_lon))

        """Collapse the time dimension to take the mean over time."""
        data_unmasked = data.collapsed('time', iris.analysis.MEAN)

        """Mask out the oceans from the data."""

        """Extract the longitudes and latitudes of the array."""
        longitude = data_unmasked.coord('longitude').points
        latitude = data_unmasked.coord('latitude').points

        """Set up Basemap to begin to mask out the ocean points from the array."""
        map = Basemap(llcrnrlon=lower_lon, llcrnrlat=lower_lat, urcrnrlon=upper_lon, urcrnrlat=upper_lat, projection='mill')
        longitude, latitude = np.meshgrid(longitude, latitude)
        x, y = map(longitude, latitude)

        """Mask out ocean points from the array, not including lakes."""

        """Create an array with the ocean data points masked."""
        masked_array = maskoceans(longitude, latitude, data_unmasked.data, resolution = 'f', grid = 1.25)

        """Convert the masked ocean data points to boolean to begin creating a boolean mask."""
        mask_bool = np.ma.filled(masked_array, 1)

        """Convert the unmasked data points to boolean to finish creating the boolean mask."""
        mask_bool[mask_bool != 1] = 0

        """Add the mask to the unmasked array."""
        data_unmasked.data = np.ma.array(data_unmasked.data, mask=mask_bool)

        """The unmasked data is now masked. Rename the array."""
        data = data_unmasked

        """Leave this code commented out. It enables printing the data on a map to test the mask to see if it has worked."""

        """
        fig = plt.figure()
        map.drawcoastlines(linewidth=1.5)
        map.drawcountries(linewidth=1)
        map.drawparallels(np.arange(-50, 60, 10), labels=[1, 0, 0, 0], fontsize=10, linewidth=0.4)
        map.drawmeridians(np.arange(-40, 80, 20), labels=[0, 0, 0, 1], fontsize=10, linewidth=0.4)
        contour_levels = np.arange(0, 110, 10)
        contour_plot = map.contourf(x, y, data_unmasked.data, contour_levels, extend='both', cmap = 'coolwarm')
        colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        fig.savefig("test1.png")
        print " plot done"
        plt.close()
        """

        """Area average the data using the iris weights method."""
        grid_areas = iris.analysis.cartography.area_weights(data)
        data = data.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

        """Print out the data for the month."""
        data = data.data
        print i
        print data

        """Append the data to the array outside the loop to produce the data for the seasonal cycle."""
        data_array = np.append(data_array, data)

    return data_array


def constrain_year(cube, time_range):
    """Constrain the years. Print the model ID, length of time dimension, and first and last model dates."""
    with iris.FUTURE.context(cell_datetime_objects=True):
        cube = cube.extract(time_range)
        time_points = cube.coord('time').points
        times = cube.coord('time').units.num2date(time_points)
        model_id = cube.attributes['model_id']
        print model_id
        print len(times)
        print times[0]
        print times[-1]
        return cube


def constrain_month(data, list_of_months):
    """Array contains the months by which the cube should be sliced."""
    time_data = num2date(data.variables['time'][:], units = data.variables['time'].units)

    list_time = np.array(time_data).tolist()

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


#seasonal_cycle(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "MIROC5", "MPI-ESM-MR", "MRI-CGCM3", "NorESM1-M/"], "amip", "evap_fraction", -10, 5, 5, 35)

#seasonal_cycle(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "MIROC5", "MPI-ESM-MR", "MRI-CGCM3", "NorESM1-M/"], "amip", "evap_fraction", -10, 5, 5, 35)
seasonal_cycle(["ACCESS1-3"], "amip", "evap_fraction", -10, 5, 5, 35)
