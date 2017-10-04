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

def seasonal_cycle(list_of_models, model_type, variable, lower_lat, upper_lat, lower_lon, upper_lon):

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

    """Find the paths to the model files themselves. """
    model_file_paths = []
    for i in directory_paths:
        files = os.listdir(i)
        for j in files:
            if "hfss" in j or "hfls" in j:
                model_file_path = os.path.join(i, j)
                model_file_paths = np.append(model_file_paths, model_file_path)
    model_file_paths.sort()
    print model_file_paths
    """Set up the figure for plotting to."""
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax2 = ax1.twinx()
    objects = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec')
    x_pos = np.arange(len(objects))
    plt.xticks(x_pos, objects)
    ax1.tick_params(axis='x', direction='in')
    ax1.tick_params(axis='y', direction='in')
    #ax2.set_ylim(0.0, 1.0)
    #ax2.tick_params(axis='y', direction='in')
    #ax1.set_zorder(ax2.get_zorder()+1)
    ax1.patch.set_visible(False)
    plt.title('Land Surface Heat Fluxes, Central Africa, AMIP, 1979-2008', fontsize = 10)

    """Define the list of colours for each model."""
    colours_array = ['red', 'darkorange', 'yellow', 'forestgreen', 'mediumblue', 'indigo', 'violet', 'black', 'darkgray']
    model_number = 0

    for j in list_of_models:

        line_colour = colours_array[model_number]
        #print line_colour
        #print model_number
        model_number +=1

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
        #print sensible_path


        if variable == 'hfsshfls':
            hfss_data_array = extract_sensible_data(sensible_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            hfls_data_array = extract_latent_data(latent_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            evap_fraction_array = []

        elif variable == 'hfss':
            hfss_data_array = extract_sensible_data(sensible_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            hfls_data_array = []
            evap_fraction_array = []

        elif variable == 'hfls':
            hfss_data_array = []
            hfls_data_array = extract_latent_data(latent_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            evap_fraction_array = []

        elif variable == 'evap_fraction':
            hfss_data_array = extract_sensible_data(sensible_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            hfls_data_array = extract_latent_data(latent_path[0], lower_lat, upper_lat, lower_lon, upper_lon)
            evap_fraction_array = hfls_data_array / (hfss_data_array + hfls_data_array)
            print evap_fraction_array

        #hfss_data_array = [48.3002731696, 52.3829183798, 52.8940286026, 46.2914752934, 43.2891164169, 45.0898602321, 47.417546342, 52.2828679393, 58.0905163153, 57.4097979651, 48.965198703, 42.840372437]
        #hfls_data_array = [70.2753424932, 71.4380859816, 80.9225648422, 88.4116364297, 84.0622804769, 74.3132192199, 66.5416929179, 65.826906716, 71.1744792892, 78.1451017016, 82.9564678829, 78.5922771205]
        #evap_fraction_array = [0.59266268, 0.57694643, 0.60472743, 0.65634443, 0.66008134, 0.62237272, 0.58390784, 0.55733665, 0.55060907, 0.57648305, 0.62883126, 0.64720878]

        #print sensible_path
        model_id = extract_model_id(sensible_path[0])

        legend = plot_data(x_pos, hfss_data_array, hfls_data_array, evap_fraction_array, ax1, line_colour, model_id, variable)

    #fig.savefig("Seasonal_Cycle_"+variable+".png")
    fig.savefig("Seasonal_Cycle_"+variable+".png", bbox_extra_artists=(legend,), bbox_inches='tight')
    print "plot done"

def plot_data(x_pos, hfss_data_array, hfls_data_array, evap_fraction_array, ax1, line_colour, model_id, variable):

    if variable == 'hfsshfls':
        ax1.plot(x_pos, hfss_data_array, zorder=1, linestyle=':', color=line_colour, label = str(model_id))
        handles, labels = ax1.get_legend_handles_labels()
        handles = [copy.copy(ha) for ha in handles]
        [ha.set_linestyle("-") for ha in handles]
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1,1))
        ax1.plot(x_pos, hfls_data_array, zorder=1, linestyle='-.', color=line_colour)
        plt.ylabel('Heat flux (W m-2)')
        plt.ylim(0, 160)
        return legend

    elif variable == 'hfss':
        ax1.plot(x_pos, hfss_data_array, zorder=1, linestyle=':', color=line_colour, label = str(model_id))
        handles, labels = ax1.get_legend_handles_labels()
        handles = [copy.copy(ha) for ha in handles]
        [ha.set_linestyle("-") for ha in handles]
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1,1))
        plt.ylabel('Sensible Heat flux (W m-2)')
        plt.ylim(0, 160)
        return legend

    elif variable == 'hfls':
        ax1.plot(x_pos, hfls_data_array, zorder=1, linestyle='-.', color=line_colour, label = str(model_id))
        handles, labels = ax1.get_legend_handles_labels()
        handles = [copy.copy(ha) for ha in handles]
        [ha.set_linestyle("-") for ha in handles]
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1,1))
        plt.ylabel('Latent Heat flux (W m-2)')
        plt.ylim(0, 160)
        return legend

    elif variable == 'evap_fraction':
        ax1.plot(x_pos, evap_fraction_array, zorder=1, linestyle='-', color=line_colour, label = str(model_id))
        handles, labels = ax1.get_legend_handles_labels()
        handles = [copy.copy(ha) for ha in handles]
        [ha.set_linestyle("-") for ha in handles]
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1,1))
        plt.ylabel('Evaporative Fraction')
        plt.ylim(0, 1)
        return legend

def extract_model_id(file_path):

    data = netCDF4.Dataset(file_path)
    model_id = data.model_id
    print model_id
    return model_id

def extract_sensible_data(file_path, lower_lat, upper_lat, lower_lon, upper_lon):

    hfss_data_array = []

    for i in np.arange(1, 13, 1):

        data = netCDF4.Dataset(file_path)

        times = extract_times(data, [i])

        hfss_data = iris.load_cube(file_path, 'surface_upward_sensible_heat_flux')
        #print hfss_data

        hfss_data = hfss_data.extract(iris.Constraint(time = times))

        hfss_data = hfss_data.intersection(latitude=(lower_lat, upper_lat), longitude=(lower_lon, upper_lon))
        #Collapse time dimension
        hfss_data_unmasked = hfss_data.collapsed('time', iris.analysis.MEAN)

        longitude = hfss_data_unmasked.coord('longitude').points

        latitude = hfss_data_unmasked.coord('latitude').points

        #print longitude
        #print latitude

        hfss = hfss_data_unmasked.data

        fig = plt.figure()

        map = Basemap(llcrnrlon=lower_lon, llcrnrlat=lower_lat, urcrnrlon=upper_lon, urcrnrlat=upper_lat, projection='mill')

        longitude, latitude = np.meshgrid(longitude, latitude)
        x, y = map(longitude, latitude)
        hfss_masked_array = maskoceans(longitude, latitude, hfss, resolution = 'f', grid = 1.25)

        mask_bool = np.ma.filled(hfss_masked_array, 1)
        # 1 is masked data, 0 is the data itself
        mask_bool[mask_bool !=1] = 0

        hfss_data = hfss_data_unmasked

        hfss_data.data = np.ma.array(hfss_data_unmasked.data, mask=mask_bool)

        #map.drawcoastlines(linewidth=1.5)
        #map.drawcountries(linewidth=1)
        #contour_levels = np.arange(0, 110, 10)
        #contour_plot = map.contourf(x, y, hfss_data.data, contour_levels, extend='both', cmap = 'coolwarm')
        #colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        #fig.savefig("test1.png")
        plt.close()

        #hfss_data.coord('latitude').guess_bounds()
        #hfss_data.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(hfss_data)
        #Take a weighted mean across lat/lon
        hfss_data = hfss_data.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

        hfss_data = hfss_data.data
        print hfss_data
        hfss_data_array = np.append(hfss_data_array, hfss_data)

    #print hfss_data_array
    return hfss_data_array

def extract_latent_data(file_path, lower_lat, upper_lat, lower_lon, upper_lon):

    hfls_data_array = []

    for i in np.arange(1, 13, 1):

        data = netCDF4.Dataset(file_path)

        times = extract_times(data, [i])

        hfls_data = iris.load_cube(file_path, 'surface_upward_latent_heat_flux')

        hfls_data = hfls_data.extract(iris.Constraint(time = times))

        hfls_data = hfls_data.intersection(latitude=(lower_lat, upper_lat), longitude=(lower_lon, upper_lon))

        hfls_data_unmasked = hfls_data.collapsed('time', iris.analysis.MEAN)

        longitude = hfls_data_unmasked.coord('longitude').points

        latitude = hfls_data_unmasked.coord('latitude').points

        #print longitude
        #print latitude

        hfls = hfls_data_unmasked.data

        fig = plt.figure()

        map = Basemap(llcrnrlon=lower_lon, llcrnrlat=lower_lat, urcrnrlon=upper_lon, urcrnrlat=upper_lat, projection='mill')

        longitude, latitude = np.meshgrid(longitude, latitude)
        x, y = map(longitude, latitude)
        hfls_masked_array = maskoceans(longitude, latitude, hfls, resolution = 'f', grid = 1.25)

        mask_bool = np.ma.filled(hfls_masked_array, 1)
        # 1 is masked data, 0 is the data itself
        mask_bool[mask_bool !=1] = 0

        hfls_data = hfls_data_unmasked

        hfls_data.data = np.ma.array(hfls_data_unmasked.data, mask=mask_bool)

        #map.drawcoastlines(linewidth=1.5)
        #map.drawcountries(linewidth=1)
        #contour_levels = np.arange(0, 110, 10)
        #contour_plot = map.contourf(x, y, hfls_data.data, contour_levels, extend='both', cmap = 'coolwarm')
        #colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        #fig.savefig("test1.png")
        plt.close()

        #hfls_data.coord('latitude').guess_bounds()
        #hfls_data.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(hfls_data)

        hfls_data = hfls_data.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

        hfls_data = hfls_data.data
        print hfls_data
        hfls_data_array = np.append(hfls_data_array, hfls_data)

    #print hfls_data_array
    return hfls_data_array

def extract_times(data, list_of_months):

    """Print the times."""

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


seasonal_cycle(["ACCESS1-3/", "bcc-csm1-1-m/", "CanAM4/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-s2/", "GFDL-CM3/", "HadGEM2-A/"], "amip", "evap_fraction", -10, 5, 5, 35)
