"""Import necessary modules for this code."""
import ensemble_cube
import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib
import matplotlib.pyplot as plt
import map_cube
from mpl_toolkits.basemap import Basemap, maskoceans
import numpy as np
import reanalysis_cube
from scipy.stats.stats import pearsonr
matplotlib.use('Agg')

def correlation(list_of_models, model_type, list_of_reanalysis, variable_x, season_name_x, variable_y, season_name_y, lower_lat, upper_lat, lower_lon, upper_lon):

    """Set up the figure for plotting to."""
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
    ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')

    """Extract the model cubes for the x variable."""
    cubes_x = map_cube.map_cube(list_of_models, model_type, variable_x, season_name_x)

    """Extract the model cubes for the y variable."""
    cubes_y = map_cube.map_cube(list_of_models, model_type, variable_y, season_name_y)

    """For each cube, mask oceans and take an area average."""
    model_variable_x = average_lat_lon(cubes_x, lower_lat, upper_lat, lower_lon, upper_lon, variable_x)
    model_variable_y = average_lat_lon(cubes_y, lower_lat, upper_lat, lower_lon, upper_lon, variable_y)

    """Define the list of colours for each model."""
    cmap = plt.get_cmap('rainbow')
    colours = [cmap(i) for i in np.linspace(0, 1, len(model_variable_x))]

    """For each model,"""
    for i in np.arange(0, len(model_variable_x), 1):

        """Select the line colour and add one to the line count variable."""
        dot_colour = colours[i]

        """Add the scatter plot. Select dot colour, marker and label."""
        ax1.scatter(model_variable_x[i], model_variable_y[i], color=dot_colour, marker='o', label=cubes_x[i].long_name)

        """Add a legend for each dot."""
        legend = plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9)

    """Extract the ensemble cube for the x variable."""
    ensemble_cube_x = ensemble_cube.ensemble(list_of_models, model_type, variable_x, season_name_x)

    """Extract the ensemble cube for the y variable."""
    ensemble_cube_y = ensemble_cube.ensemble(list_of_models, model_type, variable_y, season_name_y)

    """For each cube, mask oceans and take an area average."""
    ensemble_variable_x = average_lat_lon(iris.cube.CubeList([ensemble_cube_x]), lower_lat, upper_lat, lower_lon, upper_lon, variable_x)
    ensemble_variable_y = average_lat_lon(iris.cube.CubeList([ensemble_cube_y]), lower_lat, upper_lat, lower_lon, upper_lon, variable_x)

    """Add the scatter plot. Select dot colour, marker and label."""
    ax1.scatter(ensemble_variable_x, ensemble_variable_y, color='black', marker='o', label="Ensemble")

    """Add a legend for the ensemble mean dot."""
    legend = plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9)

    """If the x variable is precip, replace GLEAM with MSWEP."""
    if variable_x == 'pr':
        list_of_reanalysis = [i.replace("gleam", "mswep") for i in list_of_reanalysis]

    """Set up two cubes to contain the reanalysis cubes."""
    reanalysis_cubes_x_list = []
    reanalysis_cubes_y_list = []

    """For each reanalysis, extract the x cube and append to the above list."""
    for i in list_of_reanalysis:

        reanalysis_cube_x = reanalysis_cube.reanalysis([i], variable_x, season_name_x)
        reanalysis_cubes_x_list = np.append(reanalysis_cubes_x_list, reanalysis_cube_x)

    """If the y variable is precip, replace GLEAM with MSWEP."""
    if variable_y == 'pr':
        list_of_reanalysis = [i.replace("gleam", "mswep") for i in list_of_reanalysis]

    """For each reanalysis, extract the y cube and append to the above list."""
    for i in list_of_reanalysis:

        reanalysis_cube_y = reanalysis_cube.reanalysis([i], variable_y, season_name_y)
        reanalysis_cubes_y_list = np.append(reanalysis_cubes_y_list, reanalysis_cube_y)

    """For each cube, mask oceans and take an area average."""
    reanalysis_variable_x = average_lat_lon(reanalysis_cubes_x_list, lower_lat, upper_lat, lower_lon, upper_lon, variable_x)
    reanalysis_variable_y = average_lat_lon(reanalysis_cubes_y_list, lower_lat, upper_lat, lower_lon, upper_lon, variable_y)

    """Define the list of colours for each model."""
    cmap = plt.get_cmap('rainbow')
    colours = [cmap(i) for i in np.linspace(0, 1, len(reanalysis_variable_x))]

    """For each reanalysis file,"""
    for i in np.arange(0, len(reanalysis_variable_x), 1):

        """Select the line colour and add one to the line count variable."""
        cross_colour = colours[i]

        """Add the scatter plot. Select dot colour, marker and label."""
        ax1.scatter(reanalysis_variable_x[i], reanalysis_variable_y[i], color=cross_colour, marker='x', label=reanalysis_cubes_x_list[i].long_name)

        """Add a legend for each dot."""
        legend = plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9)

    """Add a line of best fit just for the models."""
    plt.plot(model_variable_x, np.poly1d(np.polyfit(model_variable_x, model_variable_y, 1))(model_variable_x), 'k')

    """Compute pearson correlation coefficient just for the models."""
    pearson = pearsonr(model_variable_x, model_variable_y)
    pearsoncoeff = pearson[0]
    pearsoncoeff_round = round(pearsoncoeff, 2)
    print str(pearsoncoeff_round)
    plt.text(0.93, 0.95, "r = "+str(pearsoncoeff_round)+"", ha='center', va='center', transform=ax1.transAxes)

    """Add labels and set x and y limits."""
    plt.xlabel('Upward Surface Latent Heat Flux (W $\mathregular{m^{-2}}$)')
    plt.ylabel('Precipitation (mm $\mathregular{day^{-1}}$)')
    plt.xlim((70, 130))
    plt.ylim((3.5, 8.0))

    """Save figure."""
    fig.savefig("Correlation.png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)


def average_lat_lon(cubelist, lower_lat, upper_lat, lower_lon, upper_lon, variable):
    """Set up a blank data array."""
    data_array = []

    """For each cube,"""
    for cube in cubelist:

        """Constrain the latitudes and longitudes of the data."""
        data_unmasked = cube.intersection(latitude=(lower_lat, upper_lat), longitude=(lower_lon, upper_lon))

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

        """The unmasked data is now masked. Rename the array as the cube."""
        cube = data_unmasked

        """Leave this code commented out. It enables printing the data on a map to test the mask to see if it has worked."""

        # fig = plt.figure()
        # map.drawcoastlines(linewidth=1.5)
        # map.drawcountries(linewidth=1)
        # map.drawparallels(np.arange(-50, 60, 10), labels=[1, 0, 0, 0], fontsize=10, linewidth=0.4)
        # map.drawmeridians(np.arange(-40, 80, 20), labels=[0, 0, 0, 1], fontsize=10, linewidth=0.4)
        # contour_levels = np.arange(0, 11, 1)
        # contour_plot = map.contourf(x, y, data_unmasked.data, contour_levels, extend='both', cmap = 'coolwarm')
        # colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        # fig.savefig("test1.png")
        # print " plot done"
        # plt.close()

        """Area average the data using the iris weights method."""
        if cube.coord('latitude').bounds is None:
            cube.coord('latitude').guess_bounds()
        if cube.coord('longitude').bounds is None:
            cube.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube)
        cube = cube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

        """Print out the data for the season."""
        data = cube.data
        print data

        """Append the data to the array outside the loop to produce the data for the correlation."""
        data_array = np.append(data_array, data)

    return data_array


correlation(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], 'amip', ["cfsr", "erai", "gleam", "jra", "merra2", "ncep-doe"], 'hfls', 'SON', 'pr', 'SON', -10, 5, 5, 35)

#correlation(["ACCESS1-3"], 'amip', ["cfsr"], 'hfls', 'SON', 'pr', 'SON', -10, 5, 5, 35)
