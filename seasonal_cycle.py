"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap, maskoceans
import netCDF4
import numpy as np
import os
import seasonal_cycle_ensemble
import seasonal_cycle_reanalysis
matplotlib.use('Agg')


def seasonal_cycle(list_of_models, model_type, variable, lower_lat, upper_lat, lower_lon, upper_lon):
    """Take the input variables, and find the paths to the relevant regridded model files."""
    """Plot these up as a seasonal cycle. Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

    if variable != "evap_fraction":

        """If variable is pr, distinguish between pr and precipitable water to find model files."""
        if variable == 'pr':
            variable = 'pr_'

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

        print model_file_paths

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
    #plt.title('Central Africa, AMIP (Land), 1979-2008', fontsize = 10)

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
        latent_path = [k for k in paths_for_this_model if 'hfls' in k]
        sensible_path = [k for k in paths_for_this_model if 'hfss' in k]
        precip_path = [k for k in paths_for_this_model if 'pr' in k]

        """Extract the model ID."""
        data = netCDF4.Dataset(paths_for_this_model[0])
        model_id = data.model_id

        """Load the cube for each variable, constrain the years and extract the seasonal cycle array."""
        """Append the seasonal cycle array to the ensemble mean array outside the loop."""

        if variable == 'pr':
            data_cube = iris.load_cube(precip_path, 'precipitation_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable)

        if variable == 'hfls':
            data_cube = iris.load_cube(latent_path, 'surface_upward_latent_heat_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable)

        if variable == 'hfss':
            data_cube = iris.load_cube(sensible_path, 'surface_upward_sensible_heat_flux')
            data_cube = constrain_year(data_cube, time_range)
            data_array = extract_seasonal_cycle_data(data_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable)

        if variable == 'evap_fraction':
            cube_latent = iris.load_cube(latent_path, 'surface_upward_latent_heat_flux')
            cube_latent = constrain_year(cube_latent, time_range)
            cube_sensible = iris.load_cube(sensible_path, 'surface_upward_sensible_heat_flux')
            cube_sensible = constrain_year(cube_sensible, time_range)

            latent_seasonal_cycle_array = extract_seasonal_cycle_data(cube_latent, latent_path[0], lower_lat, upper_lat, lower_lon, upper_lon, variable)
            sensible_seasonal_cycle_array = extract_seasonal_cycle_data(cube_sensible, sensible_path[0], lower_lat, upper_lat, lower_lon, upper_lon, variable)
            data_array = latent_seasonal_cycle_array / (latent_seasonal_cycle_array + sensible_seasonal_cycle_array, variable)

        """Add the seasonal cycle to the plot."""

        ax1.plot(x_pos, data_array, zorder=1, linestyle='-', color=line_colour, label = str(model_id))
        handles, labels = ax1.get_legend_handles_labels()
        handles[-1].set_linestyle("-")
        legend = plt.legend(handles, labels, loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9, handlelength=2.5)

        if variable == 'pr':
            plt.ylabel('Precipitation (mm $\mathregular{day^{-1}}$)')
            plt.ylim(0, 10)
        if variable == 'hfls':
            plt.ylabel('Latent Heat Flux (W $\mathregular{m^{-2}}$)')
            plt.ylim(0, 160)
        if variable == 'hfss':
            plt.ylabel('Sensible Heat Flux (W $\mathregular{m^{-2}}$)')
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

    ax1.plot(x_pos, ensemble_mean_array, zorder=2, linestyle='-', linewidth=4.0, color='black', label = "Ensemble")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("-")
    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    """Define the colours for the reanalysess."""

    colours = [cmap(i) for i in np.linspace(0, 1, 6)]

    """Add the CFSR reanalysis seasonal cycle to the plot."""

    cfsr_array = seasonal_cycle_reanalysis.seasonal_cycle_reanalysis(["cfsr"], variable, lower_lat, upper_lat, lower_lon, upper_lon)
    ax1.plot(x_pos, cfsr_array, zorder=1, linestyle='--', color = colours[0], label = "CFSR")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("--")
    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    """Add the ERAI reanalysis seasonal cycle to the plot."""

    erai_array = seasonal_cycle_reanalysis.seasonal_cycle_reanalysis(["erai"], variable, lower_lat, upper_lat, lower_lon, upper_lon)
    ax1.plot(x_pos, erai_array, zorder=1, linestyle='--', color = colours[2], label = "ERA-Interim")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("--")
    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    if variable == 'hfss':

        """Add the GLEAM reanalysis seasonal cycle to the plot."""

        gleam_array = seasonal_cycle_reanalysis.seasonal_cycle_reanalysis(["gleam"], variable, lower_lat, upper_lat, lower_lon, upper_lon)
        ax1.plot(x_pos, gleam_array, zorder=1, linestyle='--', color = colours[3], label = "GLEAM-LE")
        handles, labels = ax1.get_legend_handles_labels()
        handles[-1].set_linestyle("--")
        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    if variable == 'hfls':

        """Add the MSWEP reanalysis seasonal cycle to the plot."""

        gleam_array = seasonal_cycle_reanalysis.seasonal_cycle_reanalysis(["gleam"], variable, lower_lat, upper_lat, lower_lon, upper_lon)
        ax1.plot(x_pos, gleam_array, zorder=1, linestyle='--', color = colours[3], label = "GLEAM-LE")
        handles, labels = ax1.get_legend_handles_labels()
        handles[-1].set_linestyle("--")
        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    if variable == 'pr':

        """Add the MSWEP reanalysis seasonal cycle to the plot."""

        mswep_array = seasonal_cycle_reanalysis.seasonal_cycle_reanalysis(["mswep"], variable, lower_lat, upper_lat, lower_lon, upper_lon)
        ax1.plot(x_pos, mswep_array, zorder=1, linestyle='--', color = colours[3], label = "GLEAM (MSWEP)")
        handles, labels = ax1.get_legend_handles_labels()
        handles[-1].set_linestyle("--")
        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    """Add the JRA-55 reanalysis seasonal cycle to the plot."""

    jra_array = seasonal_cycle_reanalysis.seasonal_cycle_reanalysis(["jra"], variable, lower_lat, upper_lat, lower_lon, upper_lon)
    ax1.plot(x_pos, jra_array, zorder=1, linestyle='--', color = colours[4], label = "JRA-55")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("--")
    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    """Add the MERRA2 reanalysis seasonal cycle to the plot."""

    merra2_array = seasonal_cycle_reanalysis.seasonal_cycle_reanalysis(["merra2"], variable, lower_lat, upper_lat, lower_lon, upper_lon)
    ax1.plot(x_pos, merra2_array, zorder=1, linestyle='--', color = colours[5], label = "MERRA-2")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("--")
    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    """Add the NCEP-DOE 2 reanalysis seasonal cycle to the plot."""

    doe_array = seasonal_cycle_reanalysis.seasonal_cycle_reanalysis(["doe"], variable, lower_lat, upper_lat, lower_lon, upper_lon)
    ax1.plot(x_pos, doe_array, zorder=1, linestyle='--', color = colours[1], label = "NCEP DOE-2")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("--")
    legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

    """Save the figure."""
    fig.savefig("Seasonal_Cycle_"+variable+".png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)
    print "Plot done."


def extract_seasonal_cycle_data(input_cube, lower_lat, upper_lat, lower_lon, upper_lon, variable):
    """Extract seasonal cycle data from an iris cube given the cube, its path, and the chosen lats/lons."""
    data_array = []

    with iris.FUTURE.context(cell_datetime_objects=True):

        """For each month,"""
        for i in np.arange(1, 13, 1):

            """Constrain the data for each month only."""
            time_range = iris.Constraint(time=lambda cell: cell.point.month == i)
            data = input_cube.extract(time_range)

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
            if variable == 'pr':
                data = data*86400
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

#seasonal_cycle(["IPSL-CM5B-LR"], "amip", "pr", -10, 5, 5, 35)

seasonal_cycle(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], "amip", "pr", -10, 5, 5, 35)

#seasonal_cycle(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "MIROC5", "MPI-ESM-MR", "MRI-CGCM3", "NorESM1-M/"], "amip", "hfls", -10, 5, 5, 35)
#seasonal_cycle(["ACCESS1-3"], "amip", "hfls", -10, 5, 5, 35)
