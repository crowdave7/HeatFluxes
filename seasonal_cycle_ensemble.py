"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib
from mpl_toolkits.basemap import Basemap, maskoceans
import netCDF4
from netCDF4 import num2date
from netCDF4 import date2num
import numpy as np
import os
matplotlib.use('Agg')


def seasonal_cycle_ensemble(list_of_models, model_type, variable, lower_lat, upper_lat, lower_lon, upper_lon):
    """Take the input variables, and find the paths to the relevant regridded model files."""
    """Plot these up as a seasonal cycle. Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles"

    if variable != "evap_fraction":

        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                print path
                for j in list_of_models:
                    for char in '/':
                        j = j.replace(char,'')
                    if j in path and model_type in path and variable in path and ('hfss' in path or 'hfls' in path):
                        model_file_paths = np.append(model_file_paths, path)
        model_file_paths.sort()

    if variable == "evap_fraction":

        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                print path
                for j in list_of_models:
                    for char in '/':
                        j = j.replace(char,'')
                    if j in path and model_type in path and ('hfss' in path or 'hfls' in path):
                        model_file_paths = np.append(model_file_paths, path)
        model_file_paths.sort()

        print model_file_paths

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    """Load the data from the model file paths into a cube. Constrain the input years"""
    """Print the model ID, length of time dimension, and first and last model dates."""

    """Define a model number to begin with."""
    model_number = 0

    """Set up the array that contains the ensemble mean seasonal cycle."""
    ensemble_mean_array = np.zeros((len(list_of_models), 12))

    """For each model,"""
    for j in list_of_models:
        for char in '/':
            j = j.replace(char,'')

        """Extract the sensible and latent paths."""
        paths_for_this_model = [k for k in model_file_paths if j in k]
        latent_path = [k for k in paths_for_this_model if 'hfls' in k]
        sensible_path = [k for k in paths_for_this_model if 'hfss' in k]

        """Load the cube for each variable, constrain the years and extract the seasonal cycle array."""
        """Append the seasonal cycle array to the ensemble mean array outside the loop."""

        if variable == 'hfls':
            data_cube = iris.load_cube(latent_path, 'surface_upward_latent_heat_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, latent_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            ensemble_mean_array[model_number] = data_array

        if variable == 'hfss':
            data_cube = iris.load_cube(sensible_path, 'surface_upward_sensible_heat_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, sensible_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            ensemble_mean_array[model_number] = data_array

        if variable == 'evap_fraction':
            cube_latent = iris.load_cube(latent_path, 'surface_upward_latent_heat_flux')
            cube_latent = constrain_year(cube_latent, time_range)
            cube_sensible = iris.load_cube(sensible_path, 'surface_upward_sensible_heat_flux')
            cube_sensible = constrain_year(cube_sensible, time_range)

            latent_seasonal_cycle_array = extract_seasonal_cycle_data(cube_latent, latent_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            sensible_seasonal_cycle_array = extract_seasonal_cycle_data(cube_sensible, sensible_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            data_array = latent_seasonal_cycle_array / (latent_seasonal_cycle_array + sensible_seasonal_cycle_array)
            ensemble_mean_array[model_number] = data_array

        """Add 1 to the model number to loop through the next model."""
        model_number +=1

    """Take the mean seasonal cycle across the models."""
    ensemble_mean_array = np.mean(ensemble_mean_array, axis = 0)
    return ensemble_mean_array

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
        data.coord('latitude').guess_bounds()
        data.coord('longitude').guess_bounds()
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


#seasonal_cycle_ensemble(["CSIRO-Mk3-6-0", "HadGEM2-A"], "amip", "hfss", -10, 5, 5, 35)
