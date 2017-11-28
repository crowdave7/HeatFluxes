"""Import necessary modules for this code."""
import bowen_ratio_map_cube
import bowen_ratio_ensemble_cube
import bowen_ratio_reanalysis_cube
import ensemble_cube
import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import map_cube
from mpl_toolkits.basemap import Basemap, maskoceans, interp
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

    """If one of the input variable is evapotranspiration, use latent heat flux for now."""
    input_variable_x = variable_x
    input_variable_y = variable_y

    if input_variable_x == 'evapotranspiration':
        variable_x = 'hfls'
    if input_variable_y == 'evapotranspiration':
        variable_y = 'hfls'

    """Extract the model cubes for the x variable."""
    if input_variable_x == 'evap_fraction':
        cubes_x = bowen_ratio_map_cube.map_cube(list_of_models, model_type, season_name_x)
    else:
        cubes_x = map_cube.map_cube(list_of_models, model_type, variable_x, season_name_x)

    """Extract the model cubes for the y variable."""
    if input_variable_y == 'evap_fraction':
        cubes_y = bowen_ratio_map_cube.map_cube(list_of_models, model_type, season_name_y)
    else:
        cubes_y = map_cube.map_cube(list_of_models, model_type, variable_y, season_name_y)

    """For each cube, mask oceans and take an area average."""
    model_variable_x = average_lat_lon(cubes_x, lower_lat, upper_lat, lower_lon, upper_lon, variable_x, season_name_x)
    model_variable_y = average_lat_lon(cubes_y, lower_lat, upper_lat, lower_lon, upper_lon, variable_y, season_name_y)

    """If one of the input variables is evapotranspiration, convert latent heat flux to evapotranspiration (divide by 28.00 at 30C)"""
    if input_variable_x == 'evapotranspiration':
        model_variable_x = [float(i) / 28.00 for i in model_variable_x]

    if input_variable_y == 'evapotranspiration':
        model_variable_y = [float(i) / 28.00 for i in model_variable_y]

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
    if input_variable_x == 'evap_fraction':
        ensemble_cube_x = bowen_ratio_ensemble_cube.bowen_ensemble(list_of_models, model_type, season_name_x)
    else:
        ensemble_cube_x = ensemble_cube.ensemble(list_of_models, model_type, variable_x, season_name_x)

    """Extract the ensemble cube for the y variable."""
    if input_variable_y == 'evap_fraction':
        ensemble_cube_y = bowen_ratio_ensemble_cube.bowen_ensemble(list_of_models, model_type, season_name_y)
    else:
        ensemble_cube_y = ensemble_cube.ensemble(list_of_models, model_type, variable_y, season_name_y)

    """For each cube, mask oceans and take an area average."""
    ensemble_variable_x = average_lat_lon(iris.cube.CubeList([ensemble_cube_x]), lower_lat, upper_lat, lower_lon, upper_lon, variable_x, season_name_x)
    ensemble_variable_y = average_lat_lon(iris.cube.CubeList([ensemble_cube_y]), lower_lat, upper_lat, lower_lon, upper_lon, variable_y, season_name_y)

    """If one of the input variables is evapotranspiration, convert latent heat flux to evapotranspiration (div by 28.00 at 30C)"""
    if input_variable_x == 'evapotranspiration':
        ensemble_variable_x = [float(i) / 28.00 for i in ensemble_variable_x]

    if input_variable_y == 'evapotranspiration':
        ensemble_variable_y = [float(i) / 28.00 for i in ensemble_variable_y]

    """Add the scatter plot. Select dot colour, marker and label."""
    ax1.scatter(ensemble_variable_x, ensemble_variable_y, color='black', marker='o', label="Ensemble")

    """Add a legend for the ensemble mean dot."""
    legend = plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9)

    """Set up a cubelist for the reanalysis x variable."""
    reanalysis_cubes_x_list = []

    """If the x variable is precip, replace GLEAM with MSWEP."""
    if variable_x == 'pr':
        list_of_reanalysis = [i.replace("gleam", "mswep") for i in list_of_reanalysis]

    """For each reanalysis, extract the x cube and append to the above list."""

    for i in list_of_reanalysis:
        if input_variable_x == 'evap_fraction':
            reanalysis_cube_x = bowen_ratio_reanalysis_cube.reanalsis_bowen([i], season_name_x)
        else:
            reanalysis_cube_x = reanalysis_cube.reanalysis([i], variable_x, season_name_x)
        reanalysis_cubes_x_list = np.append(reanalysis_cubes_x_list, reanalysis_cube_x)

    print reanalysis_cubes_x_list

    """Take an average and mask the oceans."""
    reanalysis_variable_x = average_lat_lon(reanalysis_cubes_x_list, lower_lat, upper_lat, lower_lon, upper_lon, variable_x, season_name_x)

    """If variable x is precip and we have converted GLEAM to MSWEP, convert back before extracting variable y."""
    if variable_x == 'pr':
        list_of_reanalysis = [i.replace("mswep", "gleam") for i in list_of_reanalysis]

    """If the y variable is precip, replace GLEAM with MSWEP."""
    if variable_y == 'pr':
        list_of_reanalysis = [i.replace("gleam", "mswep") for i in list_of_reanalysis]

    """Set up a cubelist for the reanalysis y variable."""
    reanalysis_cubes_y_list = []

    for i in list_of_reanalysis:
        if input_variable_y == 'evap_fraction':
            reanalysis_cube_y = bowen_ratio_reanalysis_cube.reanalysis_bowen([i], season_name_y)
        else:
            reanalysis_cube_y = reanalysis_cube.reanalysis([i], variable_y, season_name_y)
        reanalysis_cubes_y_list = np.append(reanalysis_cubes_y_list, reanalysis_cube_y)

    print reanalysis_cubes_y_list

    """For each cube, mask oceans and take an area average."""
    reanalysis_variable_y = average_lat_lon(reanalysis_cubes_y_list, lower_lat, upper_lat, lower_lon, upper_lon, variable_y, season_name_y)

    """If one of the input variables is evapotranspiration, convert latent heat flux to evapotranspiration (div by 28.00 at 30C)"""
    if input_variable_x == 'evapotranspiration':
        reanalysis_variable_x = [float(i) / 28.00 for i in reanalysis_variable_x]

    if input_variable_y == 'evapotranspiration':
        reanalysis_variable_y = [float(i) / 28.00 for i in reanalysis_variable_y]

    """Define the list of colours for each model."""
    cmap = plt.get_cmap('rainbow')
    colours = [cmap(i) for i in np.linspace(0, 1, len(reanalysis_variable_x))]

    """For each reanalysis file,"""
    for i in np.arange(0, len(reanalysis_variable_x), 1):

        """Select the cross colour."""
        cross_colour = colours[i]

        """Add the scatter plot. Select cross colour, marker and label."""
        ax1.scatter(reanalysis_variable_x[i], reanalysis_variable_y[i], color=cross_colour, marker='x', label=reanalysis_cubes_x_list[i].long_name)

        """Add a legend for each dot."""
        legend = plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=9)

    # """Perform reduced major axis regression to plot the LOBF."""
    # standard_dev_variable_x = np.std(model_variable_x)
    # standard_dev_variable_y = np.std(model_variable_y)
    # gradient = standard_dev_variable_y / standard_dev_variable_x
    # mean_variable_x = np.mean(model_variable_x)
    # mean_variable_y = np.mean(model_variable_y)
    # intercept = mean_variable_y - (gradient*(mean_variable_x))
    # y_values_lobf = [gradient*i + intercept for i in model_variable_x]
    # plt.plot(model_variable_x, y_values_lobf, 'k')

    """Add a line of best fit just for the models."""
    #plt.plot(model_variable_x, np.poly1d(np.polyfit(model_variable_x, model_variable_y, 1))(model_variable_x), 'k')

    # """Compute pearson correlation coefficient just for the models."""
    # pearson = pearsonr(model_variable_x, model_variable_y)
    # pearsoncoeff = pearson[0]
    # pearsoncoeff_round = round(pearsoncoeff, 2)
    # print str(pearsoncoeff_round)
    # plt.text(0.95, 0.96, "r = "+str(pearsoncoeff_round)+"", ha='center', va='center', transform=ax1.transAxes, fontsize=8)

    """Add labels and set x and y limits."""
    if variable_x == 'hfls' and input_variable_x != 'evapotranspiration':
        plt.xlabel('Upward Surface Latent Heat Flux (W $\mathregular{m^{-2}}$)')
        plt.xlim((60, 120))
    if variable_y == 'hfls' and input_variable_y != 'evapotranspiration':
        plt.ylabel('Upward Surface Latent Heat Flux (W $\mathregular{m^{-2}}$)')
        plt.ylim((60, 130))
    if variable_x == 'pr':
        plt.xlabel('Precipitation (mm $\mathregular{day^{-1}}$)')
        plt.xlim((3.0, 8.0))
    if variable_y == 'pr':
        plt.ylabel('Precipitation (mm $\mathregular{day^{-1}}$)')
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.ylim((3.0, 8.0))
    if input_variable_x == 'evapotranspiration':
        plt.xlabel('Upward Evapotranspiration Flux (mm $\mathregular{day^{-1}}$)')
        plt.xlim((2.0, 4.5))
    if input_variable_y == 'evapotranspiration':
        plt.ylabel('Upward Evapotranspiration Flux (mm $\mathregular{day^{-1}}$)')
        plt.ylim((2.0, 4.5))
    if input_variable_x == 'mrsos':
        plt.xlabel('Volumetric Soil Moisture Content (%)')
        plt.xlim(15, 35)
    if input_variable_y == 'mrsos':
        plt.ylabel('Volumetric Soil Moisture Content (%)')
        plt.ylim(15, 35)
    if input_variable_x == 'evap_fraction':
        plt.xlabel('Evaporative Fraction')
        plt.xlim(0.5, 0.9)
    if input_variable_y == 'evap_fraction':
        plt.ylabel('Evaporative Fraction')
        plt.ylim(0.5, 0.9)
    if input_variable_x == 'tran':
        plt.xlabel('Upward Transpiration Flux (mm $\mathregular{day^{-1}}$)')
    if input_variable_y == 'tran':
        plt.ylabel('Upward Transpiration Flux (mm $\mathregular{day^{-1}}$)')
    if input_variable_x == 'evspsblsoi':
        plt.xlabel('Upward Bare Soil Evaporation Flux (mm $\mathregular{day^{-1}}$)')
    if input_variable_y == 'evspsblsoi':
        plt.ylabel('Upward Bare Soil Evaporation Flux (mm $\mathregular{day^{-1}}$)')
    if input_variable_x == 'evspsbl':
        plt.xlabel('Upward Evaporation Flux (mm $\mathregular{day^{-1}}$)')
    if input_variable_y == 'evspsbl':
        plt.ylabel('Upward Evaporation Flux (mm $\mathregular{day^{-1}}$)')
    if input_variable_x == 'evspsblveg':
        plt.ylabel('Upward Evaporation Flux from Canopy (mm $\mathregular{day^{-1}}$)')
    if input_variable_y == 'evspsblveg':
        plt.ylabel('Upward Evaporation Flux from Canopy (mm $\mathregular{day^{-1}}$)')

    """Save figure."""
    fig.savefig("Correlation_"+variable_x+"_"+variable_y+".png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)


def average_lat_lon(cubelist, lower_lat, upper_lat, lower_lon, upper_lon, variable, season_name):
    """Set up a blank data array."""
    data_array = []
    count = 0
    """For each cube,"""

    for cube in cubelist:

        """Constrain the latitudes and longitudes of the data. Set up two cubes for transposing."""
        data_unmasked = cube.intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))
        data_unmasked1 = cube.intersection(latitude=(lower_lat-2, upper_lat+1), longitude=(lower_lon-1, upper_lon+1))

        """If the coordinate needs transposing because the first coordinate is lon rather than that, transpose the data."""
        coord_names = [coord.name() for coord in data_unmasked.coords()]

        """If the first coordinate is longitude,"""
        if coord_names[0] == 'longitude':

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

        """If the first coordinate is latitude,"""
        if coord_names[0] == 'latitude':

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
        model_id = data_unmasked1.long_name
        fig.savefig("mask_"+variable+"_"+season_name+"_"+model_id+".png")
        print "plot done"

        plt.close()
        """

        """Calculate the mean of the data array (excluding nans in mask) for correlation plot."""
        data = np.nanmean(mdata)
        print data
        """Append the data to the array outside the loop to produce the data for the correlation."""
        data_array = np.append(data_array, data)

        count +=1
    return data_array

correlation(["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], 'amip', [], 'mrsos', 'SON', 'evspsblveg', 'SON', -10, 5, 5, 35)

#correlation(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "MIROC5", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], 'amip', [], 'mrsos', 'SON', 'evspsbl', 'SON', -10, 5, 5, 35)

#correlation(["bcc-csm1-1/", "BNU-ESM", "CanAM4", "GFDL-HIRAM-C360", "GISS-E2-R/", "inmcm4", "MIROC5", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], 'amip', [], 'mrsos', 'SON', 'evspsblsoi', 'SON', -10, 5, 5, 35)

#correlation(["bcc-csm1-1/", "BNU-ESM", "CanAM4", "GFDL-HIRAM-C360", "GISS-E2-R/", "inmcm4", "MIROC5", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], 'amip', [], 'mrsos', 'SON', 'tran', 'SON', -10, 5, 5, 35)

#correlation(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "MIROC5", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], 'amip', [], 'mrsos', 'SON', 'evapotranspiration', 'SON', -10, 5, 5, 35)


#correlation(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], 'amip', ["cfsr", "erai", "gleam", "jra", "merra2", "ncep-doe"], 'evaporation', 'SON', 'pr', 'SON', -10, 5, 5, 35)


#correlation(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], 'amip', ["cfsr", "erai", "gleam", "jra", "merra2", "ncep-doe"], 'hfls', 'SON', 'pr', 'SON', -10, 5, 5, 35)

#correlation(["ACCESS1-3", "bcc-csm1-1/"], 'amip', ["cfsr", "gleam"], 'hfls', 'SON', 'pr', 'SON', -10, 5, 5, 35)
