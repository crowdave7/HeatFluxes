"""Import necessary modules for this code."""
import matplotlib
import matplotlib.pyplot as plt
import netCDF4
from netCDF4 import num2date
import numpy as np
from mpl_toolkits.basemap import Basemap, maskoceans
from mpl_toolkits.basemap import shiftgrid
import os
import iris
import iris.analysis
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

    print len(model_file_paths)

    """For each model file, plot the map."""
    for i in model_file_paths:

        """Extract the data, model ID, variable name and variable units."""
        data = netCDF4.Dataset(i)
        model_id = data.model_id
        variable_name = data.variables[variable].long_name
        variable_units = data.variables[variable].units
        print model_id

        """Print a list of the variables."""
        # print data.variables

        """Find the dimensions and their sizes."""
        #for i in data.dimensions:
            #print data.dimensions[str(i)].name
            #print data.dimensions[str(i)].size

        if variable == 'hfls':
            name = 'surface_upward_latent_heat_flux'
            model_data = iris.load_cube(i, name)
        if variable == 'hfss':
            name = 'surface_upward_sensible_heat_flux'
            model_data = iris.load_cube(i, name)

        reanalysis_data = iris.load_cube("/ouce-home/data_not_backed_up/analysis/erai/1.0x1.0/daily/precip/nc/erai.totprecip.dt.1979-2013.nc")

        regridded_model_data = model_data.regrid(reanalysis_data, iris.analysis.Linear())
        #print regridded_model_data

        #Select longitudes and latitudes for Africa
        regridded_model_data = regridded_model_data.intersection(latitude=(-40, 40), longitude=(-20, 60))
        print regridded_model_data


        #Take mean for plotting, cgeck using hfls

        hfls_data_unmasked = regridded_model_data.collapsed('time', iris.analysis.MEAN)

        longitude = hfls_data_unmasked.coord('longitude').points

        latitude = hfls_data_unmasked.coord('latitude').points

        print longitude
        print latitude

        hfls = hfls_data_unmasked.data

        fig = plt.figure()

        map = Basemap(llcrnrlon=-22, llcrnrlat=-22, urcrnrlon=62, urcrnrlat=12, projection='mill')

        longitude, latitude = np.meshgrid(longitude, latitude)
        x, y = map(longitude, latitude)
        hfls_masked_array = maskoceans(longitude, latitude, hfls, resolution = 'f', grid = 1.25)

        mask_bool = np.ma.filled(hfls_masked_array, 1)
        # 1 is masked data, 0 is the data itself
        mask_bool[mask_bool !=1] = 0

        hfls_data = hfls_data_unmasked

        hfls_data.data = np.ma.array(hfls_data_unmasked.data, mask=mask_bool)

        map.drawcoastlines(linewidth=1.5)
        map.drawcountries(linewidth=1)
        contour_levels = np.arange(0, 165, 15)
        contour_plot = map.contourf(x, y, hfls_data.data, contour_levels, extend='both', cmap = 'coolwarm')
        colour_bar = map.colorbar(contour_plot, location='bottom', pad='15%')
        colour_bar.set_ticks([0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150])
        colour_bar.set_ticklabels([0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150])
        colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)
        fig.savefig("test1.png")
        plt.close()

        iris.save(regridded_model_data, "/ouce-home/students/kebl4396/Paper1/regridded/"+variable+"_"+model_id+"_regridded_Africa.nc")
        print model_id+"done"


#plot_map(["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/"], "amip", "hfls", [1,2,12], "DJF")

plot_map(["GFDL-CM3/"], "amip", "hfls", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "Climatology")
