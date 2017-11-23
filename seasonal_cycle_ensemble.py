"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans, interp
import numpy as np
import os
matplotlib.use('Agg')


def seasonal_cycle_ensemble(list_of_models, model_type, variable, lower_lat, upper_lat, lower_lon, upper_lon):
    """Take the input variables, and find the paths to the relevant regridded model files."""
    """Plot these up as a seasonal cycle. Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles"

    if variable != "evap_fraction":

        """If variable is pr, distinguish between pr and precipitable water to find model files."""
        if variable == 'pr':
            variable = 'pr_'

        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                for j in list_of_models:
                    for char in '/':
                        j = j.replace(char,'')
                    if j in path and model_type in path and variable in path:
                        model_file_paths = np.append(model_file_paths, path)

        model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

        """If variable is pr_, convert variable back to pr"""
        if variable == 'pr_':
            variable = 'pr'

    if variable == "evap_fraction":

        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                for j in list_of_models:
                    for char in '/':
                        j = j.replace(char,'')
                    if j in path and model_type in path and ('hfss' in path or 'hfls' in path):
                        model_file_paths = np.append(model_file_paths, path)

        model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

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
        precip_path = [k for k in paths_for_this_model if 'pr' in k]
        sm_path = [k for k in paths_for_this_model if 'mrsos' in k]

        """Load the cube for each variable, constrain the years and extract the seasonal cycle array."""
        """Append the seasonal cycle array to the ensemble mean array outside the loop."""

        if variable == 'pr':
            data_cube = iris.load_cube(precip_path, 'precipitation_flux')
            data_cube = constrain_year(data_cube, time_range)
            """Select the model ID"""
            model_id = data_cube.attributes['model_id']
            """Multiply precip data by 86400."""
            data_cube = iris.analysis.maths.multiply(data_cube, 86400)
            """Reassign a model ID to the new multiplied cube."""
            data_cube.long_name = model_id
            data_array = extract_seasonal_cycle_data(data_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable)
            ensemble_mean_array[model_number] = data_array

        if variable == 'hfls':
            data_cube = iris.load_cube(latent_path, 'surface_upward_latent_heat_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable)
            ensemble_mean_array[model_number] = data_array

        if variable == 'hfss':
            data_cube = iris.load_cube(sensible_path, 'surface_upward_sensible_heat_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable)
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

        if variable == 'mrsos':
            data_cube = iris.load_cube(sm_path, 'moisture_content_of_soil_layer')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable)
            ensemble_mean_array[model_number] = data_array

        """Add 1 to the model number to loop through the next model."""
        model_number +=1

    """Take the mean seasonal cycle across the models."""
    ensemble_mean_array = np.mean(ensemble_mean_array, axis = 0)
    return ensemble_mean_array

def extract_seasonal_cycle_data(input_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable):
    """Extract seasonal cycle data from an iris cube given the cube, its path, and the chosen lats/lons."""

    data_array = []

    with iris.FUTURE.context(cell_datetime_objects=True):

        """Constrain the latitudes and longitudes of the data. Set up two cubes for transposing."""
        original_data_unmasked = input_cube.intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        data_unmasked1 = input_cube.intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))

        """If the coordinate needs transposing because the first coordinate is lon rather than that, transpose the data."""
        coord_names = [coord.name() for coord in original_data_unmasked.coords()]

        print coord_names

        """For each month,"""
        for i in np.arange(1, 13, 1):
            """Constrain the data for each month only."""
            time_range = iris.Constraint(time=lambda cell: cell.point.month == i)
            data_unmasked = original_data_unmasked.extract(time_range)
            data_unmasked = data_unmasked.collapsed('time', iris.analysis.MEAN)

            """If the first coordinate is longitude,"""
            if coord_names[1] == 'longitude':

                """Set up grid of longitudes and latitudes for basemap."""
                map = Basemap(llcrnrlon=lower_lon, llcrnrlat=lower_lat, urcrnrlon=upper_lon, urcrnrlat=upper_lat, projection='mill')
                longitude = data_unmasked1.coord('longitude').points
                latitude = data_unmasked1.coord('latitude').points
                longitude, latitude = np.meshgrid(longitude, latitude)
                x, y = map(longitude, latitude)

                """Set up grid replacing each gridpoint with a 5x5 grid point."""
                x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*5)
                y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*5)
                x2, y2 = np.meshgrid(x2, y2)

                """Transpose the data to set lat first rather than lon."""
                data_unmasked = np.transpose(data_unmasked1.data, (1, 0))

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
                model_id = "Ensemble"
                fig.savefig("mask_"+variable+"_"+str(i)+"_"+model_id+".png")
                print "plot done"

                plt.close()
                """

                """Calculate the mean of the data array (excluding nans in mask) for correlation plot."""
                data = np.nanmean(mdata)

                print i
                print data

                """Append the data to the array outside the loop to produce the data for the seasonal cycle."""
                data_array = np.append(data_array, data)

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
                model_id = "Ensemble"
                fig.savefig("mask_"+variable+"_"+str(i)+"_"+model_id+".png")
                print "plot done"

                plt.close()
                """

                """Calculate the mean of the data array (excluding nans in mask) for correlation plot."""
                data = np.nanmean(mdata)

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


#seasonal_cycle_ensemble(["IPSL-CM5B-LR"], "amip", "hfss", -10, 5, 5, 35)
