#!/usr/bin/env python

"""Import necessary modules for this code."""
import cartopy as cart
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris
import iris.analysis
import iris.coord_categorisation
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import time



def build_cubelist(i, array):

    cubes_in_file = iris.load(model_file_paths[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        if len(coord_names) == 3:
            array[i] = cube
        if len(coord_names) == 4 and 'depth' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'Level' in coord_names:
            array[i] = cube

def model_file_paths_func(list_of_models, model_type, variable):

    if variable == 'nrad':
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1NetRadiationModelFiles"
    else:
        root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

    """If variable is pr, distinguish between pr and precipitable water to find model files."""
    if variable == 'pr':
        variable = 'pr_'

    """If variable is evspsbl, distinguish between evspsbl and evspsblsoi to find model files."""
    if variable == 'evspsbl':
        variable = 'evspsbl_'

    """FIND MODEL PATHS."""
    print "finding model paths"
    model_file_paths = find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory)

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """If variable is evspsbl_, convert variable back to evspsbl"""
    if variable == 'evspsbl_':
        variable = 'evspsbl'

    return model_file_paths

def find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory):
    "model_file_paths"

    if variable == 'nrad':
        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                list_of_models = [i.replace("bcc-csm1-1/", "bcc-csm1-1_") for i in list_of_models]
                list_of_models = [i.replace('/', '') for i in list_of_models]
                for j in list_of_models:
                    if j in path and model_type in path and variable in path:
                        model_file_paths = np.append(model_file_paths, path)

        model_file_paths_nrad_sorted = []

        for i in list_of_models:
            print "hi8"
            print i
            for j in model_file_paths:
                if i in j:
                    print "hi3"
                    print i
                    print j
                    print "hi4"
                    model_file_paths_nrad_sorted = np.append(model_file_paths_nrad_sorted, j)

        print model_file_paths_nrad_sorted
        model_file_paths = model_file_paths_nrad_sorted

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

    if variable not in ["nrad", "evap_fraction"]:

        """Find the paths to the directories containing the model data"""
        directory_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in directories:
                path = os.path.join(root, i)
                for j in list_of_models:
                    if j in path and model_type in path and ensemble in path:
                        print path
                        directory_paths = np.append(directory_paths, path)

        """Find the model files and their absolute paths."""
        model_file_paths = []
        for i in directory_paths:
            files = os.listdir(i)
            for j in files:
                if variable in j:
                    model_file_path = os.path.join(i, j)
                    model_file_paths = np.append(model_file_paths, model_file_path)

        #file_names = [os.path.basename(i) for i in model_file_paths]

        #model_file_paths_sorted = sorted(model_file_paths, key=lambda i: os.path.basename(i)[0])

        for i in list_of_models:
            if "CCSM4" in i:
                if variable == "evspsblveg":
                    model_file_paths = np.append(model_file_paths, "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/CCSM4/atlas_evspsblveg_Lmon_CCSM4_amip_r1i1p1_197901-201012.nc")
                if variable == "evspsblsoi":
                    model_file_paths = np.append(model_file_paths, "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/CCSM4/atlas_evspsblsoi_Lmon_CCSM4_amip_r1i1p1_197901-201012.nc")
                if variable == "tran":
                    model_file_paths = np.append(model_file_paths, "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/CCSM4/atlas_tran_Lmon_CCSM4_amip_r1i1p1_197901-201012.nc")

        model_file_paths_sorted = []

        for i in list_of_models:
            for j in model_file_paths:
                if i in j:
                    model_file_paths_sorted = np.append(model_file_paths_sorted, j)

        count = 0
        for i in model_file_paths_sorted:
            if "GISS-E2-R" in i:
                if "evspsblsoi" in i:
                    model_file_paths_sorted[count] = "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/atlas_evspsblsoi_Lmon_GISS-E2-R_amip_r1i1p1_188001-201012_correct.nc"
                    print model_file_paths_sorted[count]
            count +=1

        count = 0
        for i in model_file_paths_sorted:
            if "BNU-ESM" in i:
                if "evspsblveg" in i:
                    model_file_paths_sorted[count] = '/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles/atlas_evspsblveg_Lmon_BNU-ESM_amip_r1i1p1_197901-200812_correct.nc'
                    print model_file_paths_sorted[count]
            count +=1

        model_file_paths = model_file_paths_sorted

    print "hi2"
    print model_file_paths
    return model_file_paths

def slicing(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):

        """Slice time range."""
        if variable != "treeFrac":
            cubelist[i] = cubelist[i].extract(time_range_year)
            time_points = cubelist[i].coord('time').points
            times = cubelist[i].coord('time').units.num2date(time_points)
        model_id = cubelist[i].attributes['model_id']
        print model_id
        if variable != "treeFrac":
            print len(times)
            print times[0]
            print times[-1]

        """Slice the regridded cube down to the African domain."""
        cubelist[i] = cubelist[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

        """Select model ID."""
        if variable == 'nrad':
            model_id = cubelist[i].long_name
        else:
            model_id = cubelist[i].attributes['model_id']

        if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg':
            cubelist[i] = iris.analysis.maths.multiply(cubelist[i], 86400)

        """if variable is evapotranspiration, divide by 28"""
        if variable == 'evapotranspiration':
            cubelist[i] = iris.analysis.maths.divide(cubelist[i], 28)

        """if variable is bare soil evap, canopy evap or transpiration and model is IPSL-CM5B-LR:"""
        if variable == 'evspsblsoi' or variable == "evspsblveg" or variable == 'tran':
            if model_id == "IPSL-CM5B-LR":
                cubelist[i] = iris.analysis.maths.multiply(cubelist[i], 4)

        """Reassign model ID."""
        cubelist[i].long_name = model_id

        """If the input month is Climatology,"""
        if input_time == "Climatology":

            """Take the mean over the cube."""
            cubelist[i] = cubelist[i].collapsed('time', iris.analysis.MEAN)

        """ If the input is defined as a season,"""
        if input_time in ['DJF', 'MAM', 'JJA', 'SON']:

            if variable == 'treeFrac':
                cubelist[i] = cubelist[i].collapsed('time', iris.analysis.MEAN)

            else:
                """Create two coordinates on the cube to represent seasons and the season year."""
                iris.coord_categorisation.add_season(cubelist[i], 'time', name='clim_season')
                iris.coord_categorisation.add_season_year(cubelist[i], 'time', name='season_year')
                """Aggregate the data by season and season year."""
                seasonal_means = cubelist[i].aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)
                """Constrain the data by the input season. Decapitalise the input variable."""
                constraint = iris.Constraint(clim_season=input_time.lower())
                cube_season = seasonal_means.extract(constraint)
                """Take the mean over the cube."""
                cubelist[i] = cube_season.collapsed('time', iris.analysis.MEAN)

        """ If the input month is defined as a month,"""
        if input_time in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:

            if variable == 'treeFrac':
                cubelist[i] = cubelist[i].collapsed('time', iris.analysis.MEAN)

            else:
                """Create two coordinates on the cube to represent seasons and the season year."""
                iris.coord_categorisation.add_month(cubelist[i], 'time', name='month')
                iris.coord_categorisation.add_year(cubelist[i], 'time', name='year')

                """Aggregate the data by month."""
                monthly_means = cubelist[i].aggregated_by(['month', 'year'], iris.analysis.MEAN)

                """Constrain the data by the provided month. Decapitalise the input variable."""
                constraint = iris.Constraint(month=input_time)
                cube_monthly = monthly_means.extract(constraint)

                cubelist[i] = cube_monthly.collapsed('time', iris.analysis.MEAN)

        array[i] = cubelist[i]

def build_cubelist_ensemble(i, array):

    cubes_in_file = iris.load(model_file_paths_ensemble[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        if len(coord_names) == 3:
            array[i] = cube
        if len(coord_names) == 4 and 'depth' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'Level' in coord_names:
            array[i] = cube

def model_file_paths_ensemble_func(list_of_models, model_type, variable):

        """Import the data."""
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles"

        """If variable is pr, distinguish between pr and precipitable water to find model files."""
        if variable == 'pr':
            variable = 'pr_'

        """If variable is evspsbl, distinguish between evspsbl and evspsblsoi to find model files."""
        if variable == 'evspsbl':
            variable = 'evspsbl_'

        if variable == 'evapotranspiration':
            variable = 'hfls'

        model_file_paths = find_model_paths_ensemble(list_of_models, model_type, variable, root_directory)

        """If variable is pr_, convert variable back to pr"""
        if variable == 'pr_':
            variable = 'pr'

        """If variable is evspsbl_, convert variable back to evspsbl"""
        if variable == 'evspsbl_':
            variable = 'evspsbl'

        if variable == 'hfls':
            variable = 'evapotranspiration'

        return model_file_paths

def find_model_paths_ensemble(list_of_models, model_type, variable, root_directory):

    """Find the paths to the files containing the model data"""
    model_file_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            for j in list_of_models:
                if j == "bcc-csm1-1/":
                    j = "bcc-csm1-1_"
                for char in '/':
                    j = j.replace(char,'')
                if j in path and model_type in path and variable in path:
                    model_file_paths = np.append(model_file_paths, path)

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    return model_file_paths

def slicing_ensemble(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        if variable != "treeFrac":
            cubelist_ensemble[i] = cubelist_ensemble[i].extract(time_range_year)
            time_points = cubelist_ensemble[i].coord('time').points
            times = cubelist_ensemble[i].coord('time').units.num2date(time_points)
        model_id = cubelist_ensemble[i].attributes['model_id']
        print model_id
        if variable != "treeFrac":
            print len(times)
            print times[0]
            print times[-1]

        """Slice the regridded cube down to the African domain."""
        cubelist_ensemble[i] = cubelist_ensemble[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

        """Select model ID."""
        if variable == 'nrad':
            model_id = cubelist_ensemble[i].long_name
        else:
            model_id = cubelist_ensemble[i].attributes['model_id']

        if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg':
            cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 86400)

        """if variable is evapotranspiration, divide by 28"""
        if variable == 'evapotranspiration':
            cubelist_ensemble[i] = iris.analysis.maths.divide(cubelist_ensemble[i], 28)

        if variable == 'evspsblsoi' or variable == "evspsblveg" or variable == 'tran':
            if model_id == "IPSL-CM5B-LR":
                cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 4)

        """Reassign model ID."""
        cubelist_ensemble[i].long_name = model_id

        """If the input month is Climatology,"""
        if input_time == "Climatology":
            """Take the mean over the cube."""
            cubelist[i] = cubelist[i].collapsed('time', iris.analysis.MEAN)

        """ If the input month is defined as a season,"""
        if input_time in ['DJF', 'MAM', 'JJA', 'SON']:

            if variable == 'treeFrac':
                cubelist_ensemble[i] = cubelist_ensemble[i].collapsed('time', iris.analysis.MEAN)

            else:
                """Create two coordinates on the cube to represent seasons and the season year."""
                iris.coord_categorisation.add_season(cubelist_ensemble[i], 'time', name='clim_season')
                iris.coord_categorisation.add_season_year(cubelist_ensemble[i], 'time', name='season_year')

                """Aggregate the data by season and season year."""
                seasonal_means = cubelist_ensemble[i].aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

                """Constrain the data by the input season. Decapitalise the input variable."""
                constraint = iris.Constraint(clim_season=input_time.lower())
                cube_season = seasonal_means.extract(constraint)

                """Take the mean over the cube."""
                cubelist_ensemble[i] = cube_season.collapsed('time', iris.analysis.MEAN)

        """ If the input month is defined as a month,"""
        if input_time in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:

            if variable == 'treeFrac':
                cubelist_ensemble[i] = cubelist_ensemble[i].collapsed('time', iris.analysis.MEAN)

            else:
                """Create two coordinates on the cube to represent seasons and the season year."""
                iris.coord_categorisation.add_month(cubelist_ensemble[i], 'time', name='month')
                iris.coord_categorisation.add_year(cubelist_ensemble[i], 'time', name='year')

                """Aggregate the data by month."""
                monthly_means = cubelist_ensemble[i].aggregated_by(['month', 'year'], iris.analysis.MEAN)

                """Constrain the data by the provided month. Decapitalise the input variable."""
                constraint = iris.Constraint(month=input_time)
                cube_monthly = monthly_means.extract(constraint)

                cubelist_ensemble[i] = cube_monthly.collapsed('time', iris.analysis.MEAN)

    array[i] = cubelist_ensemble[i]


def build_cubelist_reanalysis(i, array):

    cubes_in_file = iris.load(reanalysis_file_paths[i])
    for j in np.arange(0, len(cubes_in_file)):
        cube = cubes_in_file[j]
        coord_names = [coord.name() for coord in cube.coords()]
        if len(coord_names) == 3:
            array[i] = cube

def reanalysis_file_paths_func(list_of_reanalysis, variable):
    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles"

    """If variable is pr, distinguish between pr and precipitable water to find model files."""
    if variable == 'pr':
        variable = 'pr_'

    """If variable is evspsbl, distinguish between evspsbl and evspsblsoi to find model files."""
    if variable == 'evspsbl':
        variable = 'evspsbl_'

    if variable == 'evapotranspiration':
        variable = 'hfls'

    reanalysis_file_paths = find_reanalysis_file_paths(list_of_reanalysis, variable, root_directory)

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """If variable is evspsbl_, convert variable back to evspsbl"""
    if variable == 'evspsbl_':
        variable = 'evspsbl'

    if variable == 'hfls':
        variable = 'evapotranspiration'

    return reanalysis_file_paths

def find_reanalysis_file_paths(list_of_reanalysis, variable, root_directory):

    """Find the paths to the files containing the reanalysis data"""
    reanalysis_file_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            for j in list_of_reanalysis:
                if j in path and variable in path:
                    reanalysis_file_paths = np.append(reanalysis_file_paths, path)

    reanalysis_file_paths = sorted(reanalysis_file_paths, key=lambda s: s.lower())

    return reanalysis_file_paths


def slicing_reanalysis(i, array):

    with iris.FUTURE.context(cell_datetime_objects=True):
        if variable != "treeFrac":
            cubelist_reanalysis[i] = cubelist_reanalysis[i].extract(time_range_year)
            time_points = cubelist_reanalysis[i].coord('time').points
            times = cubelist_reanalysis[i].coord('time').units.num2date(time_points)
        reanalysis_id = list_of_reanalysis[i]
        print reanalysis_id
        if variable != "treeFrac":
            print len(times)
            print times[0]
            print times[-1]

        """Slice the regridded cube down to the African domain."""
        cubelist_reanalysis[i] = cubelist_reanalysis[i].intersection(latitude=(-40, 40), longitude=(-30, 70))

        """if variable is evapotranspiration, divide by 28"""
        if variable == 'evapotranspiration':
            cubelist_reanalysis[i] = iris.analysis.maths.divide(cubelist_reanalysis[i], 28)

        """If the input month is Climatology,"""
        if input_time == "Climatology":
            """Take the mean over the cube."""
            cubelist_reanalysis[i] = cubelist_reanalysis[i].collapsed('time', iris.analysis.MEAN)

        """ If the input month is defined as a season,"""
        if input_time in ['DJF', 'MAM', 'JJA', 'SON']:

            """Create two coordinates on the cube to represent seasons and the season year."""
            iris.coord_categorisation.add_season(cubelist_reanalysis[i], 'time', name='clim_season')
            iris.coord_categorisation.add_season_year(cubelist_reanalysis[i], 'time', name='season_year')

            """Aggregate the data by season and season year."""
            seasonal_means = cubelist_reanalysis[i].aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

            """Constrain the data by the input season. Decapitalise the input variable."""
            constraint = iris.Constraint(clim_season=input_time.lower())
            reanalysis_data_season = seasonal_means.extract(constraint)

            """Take the mean over the cube."""
            cubelist_reanalysis[i] = reanalysis_data_season.collapsed('time', iris.analysis.MEAN)

        """ If the input month is defined as a month,"""
        if input_time in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:

            if variable == 'treeFrac':
                cubelist_reanalysis[i] = cubelist_reanalysis[i].collapsed('time', iris.analysis.MEAN)

            else:
                """Create two coordinates on the cube to represent seasons and the season year."""
                iris.coord_categorisation.add_month(cubelist_reanalysis[i], 'time', name='month')
                iris.coord_categorisation.add_year(cubelist_reanalysis[i], 'time', name='year')

                """Aggregate the data by month."""
                monthly_means = cubelist_reanalysis[i].aggregated_by(['month', 'year'], iris.analysis.MEAN)

                """Constrain the data by the provided month. Decapitalise the input variable."""
                constraint = iris.Constraint(month=input_time)
                cube_monthly = monthly_means.extract(constraint)

                cubelist_reanalysis[i] = cube_monthly.collapsed('time', iris.analysis.MEAN)

        """Select the reanalysis ID."""
        reanalysis_id = list_of_reanalysis[i]

        if reanalysis_id == "cfsr":
            cubelist_reanalysis[i].long_name = "CFSR"
            cubelist_reanalysis[i].rename("CFSR")
        if reanalysis_id == "erai":
            cubelist_reanalysis[i].long_name = "ERA-Interim"
            cubelist_reanalysis[i].rename("ERA-Interim")
        if reanalysis_id == "gleam":
            cubelist_reanalysis[i].long_name = "GLEAM"
            cubelist_reanalysis[i].rename("GLEAM")
        if reanalysis_id == "jra":
            cubelist_reanalysis[i].long_name = "JRA-55"
            cubelist_reanalysis[i].rename("JRA-55")
        if reanalysis_id == "merra2":
            cubelist_reanalysis[i].long_name = "MERRA-2"
            cubelist_reanalysis[i].rename("MERRA-2")
        if reanalysis_id == "mswep":
            cubelist_reanalysis[i].long_name = "MSWEP"
            cubelist_reanalysis[i].rename("MSWEP")
        if reanalysis_id == "ncep-ncar":
            cubelist_reanalysis[i].long_name = "NCEP/NCAR"
            cubelist_reanalysis[i].rename("NCEP/NCAR")
        if reanalysis_id == "ncep-doe":
            cubelist_reanalysis[i].long_name = "NCEP DOE-2"
            cubelist_reanalysis[i].rename("NCEP DOE-2")

    array[i] = cubelist_reanalysis[i]


def contour_lev_colour_map(lower_value, higher_value, interval, cmap):

    contour_levels = np.arange(lower_value, higher_value+interval, interval)
    cmap = matplotlib.cm.get_cmap(cmap)
    return contour_levels, cmap


def colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, interval):

    colour_bar.set_ticks(np.arange(lower_tick, upper_tick+interval, interval))
    colour_bar.set_ticklabels(np.arange(lower_tick, upper_tick+interval, interval))
    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0, labelsize=6)

    return colour_bar

def plot_map(cube, variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval):

    print "Setting up figure"

    """Setting up figure."""
    crs_latlon = ccrs.PlateCarree()

    """Plot the figure."""
    fig = plt.figure(figsize=(4,4))

    ax = plt.subplot(111, projection=crs_latlon)
    ax.set_extent([-22, 62, -24, 17], crs=crs_latlon)

    contour_plot = iplt.contourf(cube, contour_levels, cmap=cmap, extend='both')

    """Import coastlines and lake borders. Set the scale to 10m, 50m or 110m resolution for more detail."""
    coastline = cart.feature.NaturalEarthFeature(category='physical', name='coastline', scale='110m', facecolor='none')
    lake_borders = cart.feature.NaturalEarthFeature(category='physical', name='lakes', scale='110m', facecolor='none')

    """Import country borders."""
    shapefile = shapereader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries')
    reader = shapereader.Reader(shapefile)
    country_borders = reader.records()

    """Remove iris warning message."""
    iris.FUTURE.netcdf_promote = True

    """Plot the map using cartopy, and add map features."""
    ax.add_feature(coastline, zorder=5, edgecolor='k', linewidth=2)
    ax.add_feature(lake_borders, zorder=5, edgecolor='k', linewidth=1)
    for i in country_borders:
        ax.add_geometries(i.geometry, ccrs.PlateCarree(), edgecolor="black", facecolor="None")
    ax.add_feature(cart.feature.OCEAN, zorder=1, facecolor="w")

    """Define gridlines."""
    gridlines = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='black', linewidth=0.4, linestyle='--')
    gridlines.xlabels_top = True
    gridlines.xlabels_bottom = True
    gridlines.ylabels_left = True
    gridlines.ylabels_right = True
    gridlines.xlabel_style = {'size': 6, 'color': 'black'}
    gridlines.ylabel_style = {'size': 6, 'color': 'black'}
    gridlines.xlines = False
    gridlines.ylines = False
    gridlines.xformatter = LONGITUDE_FORMATTER
    gridlines.yformatter = LATITUDE_FORMATTER
    gridlines.xlocator = mticker.FixedLocator(np.arange(-40, 100, 20))
    gridlines.ylocator = mticker.FixedLocator(np.arange(-50, 70, 10))

    if cube.long_name == 'ACCESS1.3':
        cube.long_name = 'ACCESS1-3'

    #"""Add a title."""
    #plt.title(cube.long_name, fontsize=14)

    colourbar_axis = fig.add_axes([0.13, 0.2, 0.77, 0.02])

    """Add a colour bar, with ticks and labels."""
    colour_bar = plt.colorbar(contour_plot, colourbar_axis, orientation='horizontal')

    """Adjust ticks."""
    colour_bar = colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, tick_interval)

    if variable == 'pr':
        label = 'Precipitation (mm $\mathregular{day^{-1}}$)'
    if variable == 'hfls':
        label = 'Surface Upward Latent Heat Flux (W $\mathregular{m^{-2}}$)'
    if variable == 'evapotranspiration':
        label = 'Evapotranspiration (mm $\mathregular{day^{-1}}$)'
    if variable == 'hfss':
        label = 'Surface Upward Sensible Heat Flux (W $\mathregular{m^{-2}}$)'
    if variable == 'evap_fraction':
        label = 'Evaporative Fraction'
    if variable == 'nrad':
        label = 'Surface Net Radiation (W $\mathregular{m^{-2}}$)'
    if variable == 'mrsos':
        label = 'Soil Moisture Content of Upper Layer (mm)'
    if variable == 'mrso':
        label = 'Soil Moisture Content (mm)'
    if variable == 'tran':
        label = 'Transpiration (mm $\mathregular{day^{-1}}$)'
    if variable == 'evspsbl':
        label = 'Evaporation (mm $\mathregular{day^{-1}}$)'
    if variable == 'evspsblsoi':
        label = 'Bare Soil Evaporation (mm $\mathregular{day^{-1}}$)'
    if variable == 'evspsblveg':
        label = 'Evaporation from Canopy (mm $\mathregular{day^{-1}}$)'
    if variable == 'prveg':
        label = 'Precipitation Intercepted by Canopy (mm $\mathregular{day^{-1}}$)'
    if variable == 'mrros':
        label = 'Surface Runoff Flux (mm $\mathregular{day^{-1}}$)'
    if variable == 'lai':
        label = 'Leaf Area Index'
    if variable == 'mrro':
        label = 'Runoff Flux (mm $\mathregular{day^{-1}}$)'

    if variable == 'treeFrac':
        label = 'Tree Cover Fraction (%)'
        colour_bar.set_label(label, fontsize=6)

    else:
        print input_time
        colour_bar.set_label(input_time+" "+label, fontsize=8)

    cube_name = cube.long_name
    print cube_name

    if cube_name == "Composite":
        """Save the figure, close the plot and print an end statement."""
        print "hi1"
        if input_time in ['DJF', 'MAM', 'JJA', 'SON']:
            print season_index
            print "saving final figure"
            fig.savefig(variable+"_"+str(season_index)+"_Composite.png", dpi=600, bbox_inches='tight')
            plt.close()
            print "composite plot done"
        if input_time in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:
            print month_index
            print "saving final figure"
            fig.savefig(variable+"_"+str(month_index)+"_Composite.png", dpi=600, bbox_inches='tight')
            plt.close()
            print "composite plot done"

    if cube_name == "Composite Difference":
        """Save the figure, close the plot and print an end statement."""
        if input_time in ['DJF', 'MAM', 'JJA', 'SON']:
            print season_index
            print "saving final figure"
            fig.savefig(variable+"_"+str(season_index)+"_Composite_Difference.png", dpi=600, bbox_inches='tight')
            plt.close()
            print "composite difference plot done"
        if input_time in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:
            print month_index
            print "saving final figure"
            fig.savefig(variable+"_"+str(month_index)+"_Composite_Difference.png", dpi=600, bbox_inches='tight')
            print "composite difference plot done"
        else:
            print "saving final figure"
            fig.savefig("Composite_Difference.png", dpi=300, bbox_inches='tight')
            plt.close()

    if cube_name != "Composite" and cube_name != "Composite Difference":
        print "hi2"
        if input_time in ['DJF', 'MAM', 'JJA', 'SON']:
            print season_index
            print "saving final figure"
            fig.savefig(variable+"_"+cube.long_name+"_"+str(season_index)+".png")
            plt.close()
            print cube.long_name+" plot done"
        if input_time in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:
            print month_index
            print "saving final figure"
            fig.savefig(variable+"_"+cube.long_name+"_"+str(month_index)+".png")

def plot_subplot(cubelist, variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval, subplot_columns, subplot_rows):

    crs_latlon = ccrs.PlateCarree()

    if subplot_columns == 4 and subplot_rows == 4:
        fig = plt.figure(figsize=(8, 8))

    cube_number = 0

    for i in cubelist:
        long_name = i.long_name
        ax = plt.subplot(subplot_rows, subplot_columns, cube_number+1, projection=crs_latlon)
        ax.set_extent([-22, 62, -24, 17], crs=crs_latlon)

        contour_plot = iplt.contourf(i, contour_levels, cmap=cmap, extend='both')

        """Import coastlines and lake borders. Set the scale to 10m, 50m or 110m resolution for more detail."""
        coastline = cart.feature.NaturalEarthFeature(category='physical', name='coastline', scale='110m', facecolor='none')
        lake_borders = cart.feature.NaturalEarthFeature(category='physical', name='lakes', scale='110m', facecolor='none')

        """Import country borders."""
        shapefile = shapereader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries')
        reader = shapereader.Reader(shapefile)
        country_borders = reader.records()

        """Remove iris warning message."""
        iris.FUTURE.netcdf_promote = True

        """Plot the map using cartopy, and add map features."""
        ax.add_feature(coastline, zorder=5, edgecolor='k', linewidth=2)
        ax.add_feature(lake_borders, zorder=5, edgecolor='k', linewidth=1)
        for i in country_borders:
            ax.add_geometries(i.geometry, ccrs.PlateCarree(), edgecolor="black", facecolor="None")
        ax.add_feature(cart.feature.OCEAN, zorder=1, facecolor="w")

        """Define gridlines."""
        gridlines = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='black', linewidth=0.4, linestyle='--')
        gridlines.xlabels_top = False
        gridlines.xlabels_bottom = True
        gridlines.ylabels_left = True
        gridlines.ylabels_right = False
        gridlines.xlabel_style = {'size': 6, 'color': 'black'}
        gridlines.ylabel_style = {'size': 6, 'color': 'black'}
        gridlines.xlines = False
        gridlines.ylines = False
        gridlines.xformatter = LONGITUDE_FORMATTER
        gridlines.yformatter = LATITUDE_FORMATTER
        gridlines.xlocator = mticker.FixedLocator(np.arange(-40, 100, 20))
        gridlines.ylocator = mticker.FixedLocator(np.arange(-50, 70, 10))

        if long_name == 'ACCESS1.3':
            long_name = 'ACCESS1-3'

        """Add a title."""
        plt.title(long_name, fontsize=7)

        """Add 1 to the model number to loop through the next model."""
        cube_number +=1
        print cube_number

    if subplot_columns == 4 and subplot_rows == 4:
        #colourbar_axis = fig.add_axes([0.20, 0.07, 0.60, 0.02])
        colourbar_axis = fig.add_axes([0.21, 0.31, 0.60, 0.015])

    colour_bar = plt.colorbar(contour_plot, colourbar_axis, orientation='horizontal')

    """Adjust ticks."""
    colour_bar = colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, tick_interval)

    if variable == 'tran':
        label = 'Transpiration (mm $\mathregular{day^{-1}}$)'

    if subplot_columns == 4 and subplot_rows == 4:
        fig.subplots_adjust(left=0.05, right=0.98, bottom=0.36, top=0.96, wspace=0.4, hspace=0.55)

    print "saving final figure"
    fig.savefig("Gridded_"+variable+"_Composite.png", bbox_inches='tight')

# ------------------------------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":

    #OPTIONS:

    # JUST BASIC PLOT FOR 1 model or reanalysis.
    #composite_difference = "no"
    #composite_mean doesn't matter
    #list_of_models or #list_of_reanalysis = [1 model or 1 reanalysis]
    # time_difference = []

    # JUST BASIC PLOT BUT ONE PLOT FOR MULTIPLE MODELS
    # composite_difference = "no"
    # composite_mean = "no"
    # ensemble_maps = "no"
    #list_of_models = [more than 1 model]
    # time_difference = []

    # JUST BASIC PLOT BUT ONE PLOT FOR MULTIPLE MODELS - Ensemble
    # composite_difference = "no"
    # composite_mean = "no"
    # ensemble_maps = "yes"
    #list_of_models = [more than 1 model]
    # time_difference = []

    # COMPOSITE MEAN FOR MULTIPLE MODELS
    # composite_difference = "no"
    # composite_mean = "yes"
    #list_of_models = [more than 1 model]
    # time_difference = []

    # COMPOSITE DIFFERENCE FOR MULTIPLE MODELS
    # composite_difference = "yes"
    # doesn't matter about composite mean
    # fill in #list_of_models_higher_composite and #list_of_models_lower_composite
    # time_difference = []

    # COMPOSITE DIFFERENCE E.G. NOV - SEP ALL MODELS
    # composite_difference = "yes"
    # composite_mean = "yes"
    #time_difference = ["Nov", "Mar"]

    composite_difference = "yes"
    composite_mean = "yes"
    #time_diff = (2 - 1)
    time_difference = ['Mar', 'Nov']
    #time_difference = []
    label_for_time_difference = str('Nov - Mar')
    ensemble_maps = "yes"
    multiple_models = 'yes'

    # CANNOT PRINT COMPOSITE MAPS OF BOTH MODEL AND REANALYSIS AT ONCE.
    # CANNOT HAVE MORE THAN 1 REANALYSIS.
    #list_of_models = ["MRI-AGCM3-2H/"]
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/",  "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["GISS-E2-R/", "inmcm4/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_models = ['bcc-csm1-1/', 'IPSL-CM5B-LR/', 'FGOALS-g2/', 'CMCC-CM/', 'MRI-CGCM3/']
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    model_type = "amip"
    list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ["cfsr"]
    #list_of_reanalysis = ["cfsr", "erai", "gleam", "jra", "merra2", "ncep-doe"]

    #list_of_models_higher_composite = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models_lower_composite = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]


    #list_of_times = ['SON']
    # ONLY IF TIME DIFFERENCE = []
    list_of_times = ['Nov']
    variable = "tran"

    """
    djf strongest
    'inmcm4/', 'BNU-ESM/', 'CCSM4/', 'ACCESS1-3/', 'GISS-E2-R/'

    djf weakest
    'bcc-csm1-1/', 'CSIRO-Mk3-6-0/', 'CMCC-CM/', 'CNRM-CM5/', 'FGOALS-g2/' redo + dif

    mam strongest
    'inmcm4/', 'MIROC5/', 'GFDL-CM3/', 'ACCESS1-0/', 'ACCESS1-3/'

    mam weakest
    'EC-EARTH/', 'CNRM-CM5/', 'MRI-CGCM3/', 'bcc-csm1-1/', 'FGOALS-g2/' redo + dif

    jja strongest
    'MIROC5/', 'BNU-ESM/', 'ACCESS1-3/', 'ACCESS1-0/', 'CCSM4/'

    jja weakest
    'EC-EARTH/', 'bcc-csm1-1-m/', 'bcc-csm1-1/','FGOALS-g2/', 'FGOALS-s2/' no need to redo

    son strongest
    'inmcm4/', 'MIROC5/', 'BNU-ESM/', 'ACCESS1-3/', 'GFDL-HIRAM-C360/'

    son weakest
    'bcc-csm1-1/', 'IPSL-CM5B-LR/', 'FGOALS-g2/', 'CMCC-CM/', 'MRI-CGCM3/' redo + dif
    """
    lower_year = 1979
    upper_year = 2008
    lower_value = 0
    higher_value = 100
    value_interval = 10
    lower_tick = 0
    upper_tick = 100
    tick_interval = 10

    # lower_value_diff = -2.0
    # higher_value_diff = 2.25
    # value_interval_diff = 0.25
    # lower_tick_diff = -2.0
    # upper_tick_diff = 2.25
    # tick_interval_diff = 0.5

    lower_value_diff = -1.0
    higher_value_diff = 1.0
    value_interval_diff = 0.1
    lower_tick_diff = -1.0
    upper_tick_diff = 1.0
    tick_interval_diff = 0.1

    subplot_columns = 4
    subplot_rows = 4

    time_range_year = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)

    if composite_difference == "yes":
        cmap = "RdBu_r"
    else:
        cmap = "YlGnBu"

# ------------------------------------------------------------------------------------------------------------------------------------------

    if len(time_difference) == 0:

        for i in np.arange(0, len(list_of_times)):

            input_time = list_of_times[i]

            if input_time in ['DJF', 'MAM', 'JJA', 'SON']:

                if input_time == "DJF":
                    season_index = "S1"
                if input_time == "MAM":
                    season_index = "S2"
                if input_time == "JJA":
                    season_index = "S3"
                if input_time == "SON":
                    season_index = "S4"

                print input_time

            if input_time in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:

                month_index = i+1

                print input_time

                print month_index

    # ------------------------------------------------------------------------------------------------------------------------------------------

            """IF WE ARE NOT DOING A COMPOSITE DIFFERENCE (i.e. just a composite plot)"""
            if composite_difference != "yes":

        # ------------------------------------------------------------------------------------------------------------------------------------------

                """IF WE JUST WANT 1 MODEL PLOTTED."""
                if len(list_of_models) == 0:
                    pass

                if len(list_of_models) == 1:
                    """Extract the model file paths."""
                    start_time = time.time()
                    model_file_paths = model_file_paths_func(list_of_models, model_type, variable)
                    print time.time() - start_time, "seconds"

                    """Build a list of cubes from the model file paths."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(model_file_paths)):
                        p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist = array.values()

                    """Define original long name and original units for later."""
                    original_long_name = cubelist[0].long_name
                    original_unit = cubelist[0].units

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist)):
                        p = multiprocessing.Process(target=slicing, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist = array.values()

                    model_cube = cubelist[0]

                    contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
                    cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

                    plot_map(model_cube, variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval)

        # ------------------------------------------------------------------------------------------------------------------------------------------
                """IF WE JUST WANT 1 REANALYSIS PLOTTED."""

                if len(list_of_reanalysis) == 0:
                    pass

                if len(list_of_reanalysis) == 1:

                    """Build a list of cubes from the reanalysis file paths."""
                    reanalysis_file_paths = reanalysis_file_paths(list_of_reanalysis, variable)
                    """Load the cubes."""

                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(reanalysis_file_paths)):
                        p = multiprocessing.Process(target=build_cubelist_reanalysis, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_reanalysis = array.values()

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist_reanalysis)):
                        p = multiprocessing.Process(target=slicing_reanalysis, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_reanalysis = array.values()

                    reanalysis_cube = cubelist_reanalysis[0]

                    contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
                    cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

                    plot_map(reanalysis_cube, variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval)


        # ------------------------------------------------------------------------------------------------------------------------------------------

                """IF WE ARE PLOTTING SEVERAL MAPS FOR MULTIPLE MODELS."""

                if len(list_of_models) > 1 and composite_mean == "no" and ensemble_maps != 'yes':

                    """Extract the model file paths."""
                    start_time = time.time()
                    model_file_paths = model_file_paths_func(list_of_models, model_type, variable)
                    print time.time() - start_time, "seconds"

                    """Build a list of cubes from the model file paths."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(model_file_paths)):
                        p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist = array.values()

                    """Define original long name and original units for later."""
                    original_long_name = cubelist[0].long_name
                    original_unit = cubelist[0].units

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist)):
                        p = multiprocessing.Process(target=slicing, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist = array.values()

                    print cubelist

                    contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
                    cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

                    """For each model, plot."""

                    for i in np.arange(0, len(cubelist)):

                        plot_map(cubelist[i], variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval)

        # ------------------------------------------------------------------------------------------------------------------------------------------

                """IF WE ARE PLOTTING SEVERAL MAPS FOR MULTIPLE MODELS - MAPS FOR THE ENSEMBLE MEAN."""

                if len(list_of_models) > 1 and composite_mean == "no" and ensemble_maps == 'yes':

                    """Extract the model file paths."""
                    start_time = time.time()
                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)
                    print time.time() - start_time, "seconds"

                    """Build a list of cubes from the regridded model file paths."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(model_file_paths_ensemble)):
                        p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_ensemble = array.values()

                    """Define original long name and original units for later."""
                    original_long_name = cubelist_ensemble[0].long_name
                    original_unit = cubelist_ensemble[0].units

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist_ensemble)):
                        p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_ensemble = array.values()

                    print cubelist_ensemble

                    contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
                    cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

                    """For each model, plot."""

                    for i in np.arange(0, len(cubelist_ensemble)):

                        plot_map(cubelist_ensemble[i], variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval)

        # ------------------------------------------------------------------------------------------------------------------------------------------

                """ IF WE ARE PLOTTING SEVERAL MAPS FOR SEVERAL REANALYSES."""

                if len(list_of_reanalysis) > 1 and composite_mean == "no":

                    """Build a list of cubes from the reanalysis file paths."""
                    reanalysis_file_paths = reanalysis_file_paths(list_of_reanalysis, variable)
                    """Load the cubes."""

                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(reanalysis_file_paths)):
                        p = multiprocessing.Process(target=build_cubelist_reanalysis, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_reanalysis = array.values()

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist_reanalysis)):
                        p = multiprocessing.Process(target=slicing_reanalysis, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_reanalysis = array.values()

                    contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
                    cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

                    for i in np.arange(0, len(cubelist_reanalysis)):

                        plot_map(cubelist_reanalysis[i], variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval)


        # ------------------------------------------------------------------------------------------------------------------------------------------

                """IF WE ARE DOING A COMPOSITE MEAN FOR MULTIPLE MODELS"""

                if len(list_of_models) > 1 and composite_mean == "yes":

                    """Extract the regridded model file paths for the ensemble mean."""
                    print list_of_models
                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

                    """Build a list of cubes from the regridded model file paths."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(model_file_paths_ensemble)):
                        p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_ensemble = array.values()

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist_ensemble)):
                        p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_ensemble = array.values()

                    cube_1 = cubelist_ensemble[0]

                    for i in range(len(cubelist_ensemble)):
                        array = np.ma.filled(cubelist_ensemble[i].data)
                        array[array==1e20] = np.nan
                        cubelist_ensemble[i].data = array
                        cubelist_ensemble[i] = cubelist_ensemble[i].data

                    ensemble_cube = np.nansum(cubelist_ensemble, axis = 0) / float(len(cubelist_ensemble))



                    # # cube_1.data = cube_data
                    # # cube_1.long_name = "Composite"
                    # # ensemble_mean_cube = cube_1
                    # #
                    # # print ensemble_mean_cube
                    #
                    # # qplt.contourf(ensemble_mean_cube, 25, cmap="YlGnBu")
                    # # plt.gca().coastlines()
                    # # plt.show()
                    #
                    # contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
                    # cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]
                    #
                    # plot_map(ensemble_mean_cube, variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval)

        # ------------------------------------------------------------------------------------------------------------------------------------------

            """IF WE ARE DOING A COMPOSITE DIFFERENCE BETWEEN COMPOSITE MEANS"""

            if composite_difference == "yes":

                if len(list_of_models_higher_composite) == 0:
                    pass

                if len(list_of_models_higher_composite) == 1:

                    """Extract the model file paths."""
                    start_time = time.time()
                    model_file_paths = model_file_paths(list_of_models_higher_composite, model_type, variable)
                    print time.time() - start_time, "seconds"

                    """Build a list of cubes from the model file paths."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(model_file_paths)):
                        p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist = array.values()

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist)):
                        p = multiprocessing.Process(target=slicing, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist = array.values()

                    higher_composite_cube = cubelist[0]

                if len(list_of_models_higher_composite) > 1:

                    """Extract the regridded model file paths for the ensemble mean."""
                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models_higher_composite, model_type, variable)

                    print model_file_paths_ensemble

                    """Build a list of cubes from the regridded model file paths."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(model_file_paths_ensemble)):
                        p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_ensemble = array.values()

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist_ensemble)):
                        p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_ensemble = array.values()

                    """Compute the ensemble mean cube."""
                    higher_composite_cube = sum(cubelist_ensemble) / float(len(cubelist_ensemble))
                    higher_composite_cube.long_name = "Strong Composite"

                    print higher_composite_cube

                if len(list_of_models_lower_composite) == 0:
                    pass

                if len(list_of_models_lower_composite) == 1:
                    """Extract the model file paths."""

                    start_time = time.time()
                    model_file_paths = model_file_paths(list_of_models_lower_composite, model_type, variable)
                    print time.time() - start_time, "seconds"

                    """Build a list of cubes from the model file paths."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(model_file_paths)):
                        p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist = array.values()

                    """Define original long name and original units for later."""
                    original_long_name = cubelist[0].long_name
                    original_unit = cubelist[0].units

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist)):
                        p = multiprocessing.Process(target=slicing, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist = array.values()

                    lower_composite_cube = cubelist[0]

                if len(list_of_models_lower_composite) > 1:

                    print list_of_models
                    print list_of_models_lower_composite

                    """Extract the regridded model file paths for the ensemble mean."""
                    model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models_lower_composite, model_type, variable)

                    """Build a list of cubes from the regridded model file paths."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(model_file_paths_ensemble)):
                        p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_ensemble = array.values()

                    """Slice each cube by the time range."""
                    manager = multiprocessing.Manager()
                    array = manager.dict()
                    jobs = []
                    for i in np.arange(0, len(cubelist_ensemble)):
                        p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                        jobs.append(p)
                        p.start()
                    for process in jobs:
                        process.join()
                    cubelist_ensemble = array.values()

                    """Compute the ensemble mean cube."""
                    lower_composite_cube = sum(cubelist_ensemble) / float(len(cubelist_ensemble))
                    lower_composite_cube.long_name = "Weak Composite"

                    print lower_composite_cube

                qplt.contourf(higher_composite_cube, 25, cmap="YlGnBu")
                plt.gca().coastlines()
                plt.show()

                qplt.contourf(lower_composite_cube, 25, cmap="YlGnBu")
                plt.gca().coastlines()
                plt.show()

                composite_difference_cube = higher_composite_cube - lower_composite_cube
                composite_difference_cube.long_name = "Composite Difference"

                qplt.contourf(composite_difference_cube, 25, cmap="RdBu_r")
                plt.gca().coastlines()
                plt.show()

                contour_levels = contour_lev_colour_map(lower_value_diff, higher_value_diff, value_interval_diff, cmap)[0]
                cmap = contour_lev_colour_map(lower_value_diff, higher_value_diff, value_interval_diff, cmap)[1]

                plot_map(composite_difference_cube, variable, input_time, contour_levels, cmap, lower_tick_diff, upper_tick_diff, tick_interval_diff)


    #---------------------------------------------------------------------------------------------------

    """IF WE ARE DOING A COMPOSITE DIFFERENCE BETWEEN COMPOSITE MEANS OF SAME MODELS"""
    """E.g. Nov transpiration - Mar transpiration, same models"""

    if len(time_difference) != 0:

        #---------------------------------------------------------------------------------------------

        """FIRST CUBE (e.g. Nov transpiration, averaged across all models"""

        if len(list_of_models) >= 1 and composite_mean == "yes" and multiple_models == 'no':

            input_time = time_difference[0]

            """Extract the regridded model file paths for the ensemble mean."""
            print list_of_models
            model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

            """Build a list of cubes from the regridded model file paths."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(model_file_paths_ensemble)):
                p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            cube_1 = cubelist_ensemble[0]

            for i in range(len(cubelist_ensemble)):

                array = np.ma.filled(cubelist_ensemble[i].data)

                array[array==1e20] = np.nan

                cubelist_ensemble[i].data = array

                cubelist_ensemble[i] = cubelist_ensemble[i].data

            cube_data = np.nansum(cubelist_ensemble, axis = 0) / float(len(cubelist_ensemble))

            cube_1.data = cube_data

            cube_1.long_name = "Cube 1"

            print cube_1

            #---------------------------------------------------------------------------------------------

            """SECOND CUBE (e.g. Transpiration Mar, same models)"""

            input_time = time_difference[1]

            """Extract the regridded model file paths for the ensemble mean."""
            print list_of_models
            model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

            """Build a list of cubes from the regridded model file paths."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(model_file_paths_ensemble)):
                p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            """Slice each cube by the time range."""
            manager = multiprocessing.Manager()
            array = manager.dict()
            jobs = []
            for i in np.arange(0, len(cubelist_ensemble)):
                p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_ensemble = array.values()

            cube_2 = cubelist_ensemble[0]

            for i in range(len(cubelist_ensemble)):

                array = np.ma.filled(cubelist_ensemble[i].data)

                array[array==1e20] = np.nan

                cubelist_ensemble[i].data = array

                cubelist_ensemble[i] = cubelist_ensemble[i].data

            cube_data = np.nansum(cubelist_ensemble, axis = 0) / float(len(cubelist_ensemble))

            cube_2.data = cube_data

            cube_2.long_name = "Cube 2"

            composite_difference_cube = cube_1 - cube_2
            composite_difference_cube.long_name = "Composite Difference"

            print cube_1
            print cube_2
            print composite_difference_cube

            # qplt.contourf(composite_difference_cube, 25)
            #
            # plt.gca().coastlines()
            #
            # plt.show()

            contour_levels = contour_lev_colour_map(lower_value_diff, higher_value_diff, value_interval_diff, cmap)[0]
            cmap = contour_lev_colour_map(lower_value_diff, higher_value_diff, value_interval_diff, cmap)[1]

            input_time = label_for_time_difference

            plot_map(composite_difference_cube, variable, input_time, contour_levels, cmap, lower_tick_diff, upper_tick_diff, tick_interval_diff)


    #---------------------------------------------------------------------------------------------------

    """COMPOSITE DIFFERENCE, NOV-MAR TRAN, MULTIPLE MODELS, SUBPLOTTED, WITH ENSEMBLE MEAN"""

    if len(time_difference) != 0:

        if len(list_of_models) >= 1 and composite_mean == "yes" and multiple_models == 'yes':

            #--------------------------------------------------------------------------------------------

            """FIRST MONTH."""

            input_time = time_difference[0]
            if len(list_of_models) == 0:
                pass
            else:

                """Extract the model file paths."""
                start_time = time.time()
                model_file_paths = model_file_paths_func(list_of_models, model_type, variable)
                print time.time() - start_time, "seconds"

                """Build a list of cubes from the model file paths."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(model_file_paths)):
                    p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist = array.values()

                """Slice each cube by the time range."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(cubelist)):
                    p = multiprocessing.Process(target=slicing, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist = array.values()

            """ENSEMBLE MEAN."""

            if ensemble_maps == "yes":

                """Extract the regridded model file paths for the ensemble mean."""
                model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

                """Build a list of cubes from the regridded model file paths."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(model_file_paths_ensemble)):
                    p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_ensemble = array.values()

                """Slice each cube by year, month and spatial dimension."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(cubelist_ensemble)):
                    p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_ensemble = array.values()

                """Compute the ensemble mean cube, ignoring fill values."""

                cube_1 = cubelist_ensemble[0]

                for i in range(len(cubelist_ensemble)):
                    array = np.ma.filled(cubelist_ensemble[i].data)
                    array[array==1e20] = np.nan
                    cubelist_ensemble[i].data = array
                    cubelist_ensemble[i] = cubelist_ensemble[i].data

                cube_data = np.nansum(cubelist_ensemble, axis = 0) / float(len(cubelist_ensemble))

                cube_1.data = cube_data
                cube_1.long_name = "Ensemble"
                ensemble_mean_cube = cube_1

                """Append ensemble mean to the cubelist."""
                cubelist = np.append(cubelist, ensemble_mean_cube)

            if len(list_of_reanalysis) == 0:
                pass

            else:

                list_of_reanalysis.sort()
                """Build a list of cubes from the reanalysis file paths."""
                reanalysis_file_paths = reanalysis_file_paths_func(list_of_reanalysis, variable)

                """Load the cubes."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(reanalysis_file_paths)):
                    p = multiprocessing.Process(target=build_cubelist_reanalysis, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_reanalysis = array.values()

                """Slice each cube by year, month and spatial dimension."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                print len(cubelist_reanalysis)
                print np.arange(0, len(cubelist_reanalysis))
                for i in np.arange(0, len(cubelist_reanalysis)):
                    p = multiprocessing.Process(target=slicing_reanalysis, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_reanalysis = array.values()

                """Append reanalysis cubes to the cubelist."""
                cubelist = np.append(cubelist, cubelist_reanalysis)

            cubelist_1 = cubelist

            #--------------------------------------------------------------------------------------------

            """SECOND MONTH"""

            input_time = time_difference[1]

            if len(list_of_models) == 0:
                pass
            else:

                """Extract the model file paths."""
                start_time = time.time()
                model_file_paths = model_file_paths_func(list_of_models, model_type, variable)
                print time.time() - start_time, "seconds"

                """Build a list of cubes from the model file paths."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(model_file_paths)):
                    p = multiprocessing.Process(target=build_cubelist, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist = array.values()

                """Slice each cube by the time range."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(cubelist)):
                    p = multiprocessing.Process(target=slicing, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist = array.values()

            """ENSEMBLE MEAN."""

            if ensemble_maps == "yes":

                """Extract the regridded model file paths for the ensemble mean."""
                model_file_paths_ensemble = model_file_paths_ensemble_func(list_of_models, model_type, variable)

                """Build a list of cubes from the regridded model file paths."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(model_file_paths_ensemble)):
                    p = multiprocessing.Process(target=build_cubelist_ensemble, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_ensemble = array.values()

                """Slice each cube by year, month and spatial dimension."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(cubelist_ensemble)):
                    p = multiprocessing.Process(target=slicing_ensemble, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_ensemble = array.values()

                """Compute the ensemble mean cube, ignoring fill values."""

                cube_1 = cubelist_ensemble[0]

                for i in range(len(cubelist_ensemble)):
                    array = np.ma.filled(cubelist_ensemble[i].data)
                    array[array==1e20] = np.nan
                    cubelist_ensemble[i].data = array
                    cubelist_ensemble[i] = cubelist_ensemble[i].data

                cube_data = np.nansum(cubelist_ensemble, axis = 0) / float(len(cubelist_ensemble))

                cube_1.data = cube_data
                cube_1.long_name = "Ensemble"
                ensemble_mean_cube = cube_1

                """Append ensemble mean to the cubelist."""
                cubelist = np.append(cubelist, ensemble_mean_cube)

            if len(list_of_reanalysis) == 0:
                pass

            else:

                list_of_reanalysis.sort()
                """Build a list of cubes from the reanalysis file paths."""
                reanalysis_file_paths = reanalysis_file_paths_func(list_of_reanalysis, variable)

                """Load the cubes."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                for i in np.arange(0, len(reanalysis_file_paths)):
                    p = multiprocessing.Process(target=build_cubelist_reanalysis, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_reanalysis = array.values()

                """Slice each cube by year, month and spatial dimension."""
                manager = multiprocessing.Manager()
                array = manager.dict()
                jobs = []
                print cubelist_reanalysis
                print len(cubelist_reanalysis)
                print np.arange(0, len(cubelist_reanalysis))
                for i in np.arange(0, len(cubelist_reanalysis)):
                    p = multiprocessing.Process(target=slicing_reanalysis, args=(i, array))
                    jobs.append(p)
                    p.start()
                for process in jobs:
                    process.join()
                cubelist_reanalysis = array.values()

                """Append reanalysis cubes to the cubelist."""
                cubelist = np.append(cubelist, cubelist_reanalysis)

            cubelist_2 = cubelist

            print cubelist_1
            print cubelist_2

            """Collect list of model names."""

            list_of_cube_names = []

            for i in np.arange(0, len(cubelist_2)):

                cube_name = cubelist_2[i].long_name

                print cube_name

                list_of_cube_names = np.append(list_of_cube_names, cube_name)

            print list_of_cube_names

            """Subtract cubelists."""

            composite_cubelist = cubelist_2 - cubelist_1

            """Reassign model names."""

            for i in np.arange(0, len(composite_cubelist)):

                long_name = list_of_cube_names[i]

                composite_cubelist[i].long_name = long_name

            print composite_cubelist

            contour_levels = contour_lev_colour_map(lower_value_diff, higher_value_diff, value_interval_diff, cmap)[0]
            cmap = contour_lev_colour_map(lower_value_diff, higher_value_diff, value_interval_diff, cmap)[1]

            input_time = label_for_time_difference

            plot_subplot(composite_cubelist, variable, input_time, contour_levels, cmap, lower_tick_diff, upper_tick_diff, tick_interval_diff, subplot_columns, subplot_rows)
