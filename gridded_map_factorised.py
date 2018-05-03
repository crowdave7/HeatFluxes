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

    """If variable is mrso, distinguish between mrso and mrsos to find model files."""
    if variable == 'mrso':
        variable = 'mrso_'

    """If variable is mrro, distinguish between mrro and mrros to find model files."""
    if variable == 'mrro':
        variable = 'mrro_'

    if variable == 'evapotranspiration':
        variable = 'hfls'

    """FIND MODEL PATHS."""
    print "finding model paths"
    model_file_paths = find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory)

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """If variable is evspsbl_, convert variable back to evspsbl"""
    if variable == 'evspsbl_':
        variable = 'evspsbl'

    """If variable is mrso, convert variable back to mrso"""
    if variable == 'mrso_':
        variable = 'mrso'

    """If variable is mrro, convert variable back to mrro"""
    if variable == 'mrro_':
        variable = 'mrro'

    if variable == 'hfls':
        variable = 'evapotranspiration'

    return model_file_paths

def find_model_file_paths(list_of_models, model_type, ensemble, variable, root_directory):
    "model_file_paths"

    if variable == 'nrad':
        """Find the paths to the files containing the model data"""
        model_file_paths = []
        for root, directories, files in os.walk(root_directory):
            for i in files:
                path = os.path.join(root, i)
                print path
                for j in list_of_models:
                    if j == "bcc-csm1-1/":
                        j = "bcc-csm1-1_"
                    for char in '/':
                        j = j.replace(char,'')
                    if j in path and model_type in path and variable in path:
                        model_file_paths = np.append(model_file_paths, path)

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

        if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg' or variable == 'mrro' or variable == 'mrros':
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

        """ If the input month is defined as the whole year,"""
        if input_time == 'Climatology':

            """Take the mean over the cube."""
            cubelist[i] = cubelist[i].collapsed('time', iris.analysis.MEAN)

        """ If the input month is defined as a season,"""
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

        """If variable is mrso, distinguish between mrso and mrsos to find model files."""
        if variable == 'mrso':
            variable = 'mrso_'

        """If variable is mrro, distinguish between mrro and mrros to find model files."""
        if variable == 'mrro':
            variable = 'mrro_'

        if variable == 'evapotranspiration':
            variable = 'hfls'

        model_file_paths = find_model_paths_ensemble(list_of_models, model_type, variable, root_directory)

        """If variable is pr_, convert variable back to pr"""
        if variable == 'pr_':
            variable = 'pr'

        """If variable is evspsbl_, convert variable back to evspsbl"""
        if variable == 'evspsbl_':
            variable = 'evspsbl'

        """If variable is mrso, convert variable back to mrso"""
        if variable == 'mrso_':
            variable = 'mrso'

        """If variable is mrro, convert variable back to mrro"""
        if variable == 'mrro_':
            variable = 'mrro'

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

        if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg' or variable == 'mrro' or variable == 'mrros' or variable == 'mrso':
            cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 86400)

        """if variable is evapotranspiration, divide by 28"""
        if variable == 'evapotranspiration':
            cubelist_ensemble[i] = iris.analysis.maths.divide(cubelist_ensemble[i], 28)

        if variable == 'evspsblsoi' or variable == "evspsblveg" or variable == 'tran':
            if model_id == "IPSL-CM5B-LR":
                cubelist_ensemble[i] = iris.analysis.maths.multiply(cubelist_ensemble[i], 4)

        """Reassign model ID."""
        cubelist_ensemble[i].long_name = model_id

        """ If the input month is defined as the whole year,"""
        if input_time == 'Climatology':

            """Take the mean over the cube."""
            cubelist_ensemble[i] = cubelist_ensemble[i].collapsed('time', iris.analysis.MEAN)

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
        if len(coord_names) == 4 and 'depth' in coord_names:
            array[i] = cube
        if len(coord_names) == 4 and 'Level' in coord_names:
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

    print reanalysis_file_paths
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

        """ If the input month is defined as the whole year,"""
        if input_time == 'Climatology':

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
            cubelist_reanalysis[i].long_name = "GLEAM-LE"
            cubelist_reanalysis[i].rename("GLEAM-LE")
        if reanalysis_id == "jra":
            cubelist_reanalysis[i].long_name = "JRA-55"
            cubelist_reanalysis[i].rename("JRA-55")
        if reanalysis_id == "merra2":
            cubelist_reanalysis[i].long_name = "MERRA-2"
            cubelist_reanalysis[i].rename("MERRA-2")
        if reanalysis_id == "mswep":
            cubelist_reanalysis[i].long_name = "GLEAM (MSWEP)"
            cubelist_reanalysis[i].rename("GLEAM (MSWEP)")
        if reanalysis_id == "ncep-ncar":
            cubelist_reanalysis[i].long_name = "NCEP/NCAR"
            cubelist_reanalysis[i].rename("NCEP/NCAR")
        if reanalysis_id == "ncep-doe":
            cubelist_reanalysis[i].long_name = "NCEP DOE-2"
            cubelist_reanalysis[i].rename("NCEP DOE-2")

    array[i] = cubelist_reanalysis[i]


def contour_lev_colour_map(lower_value, higher_value, interval, cmap):

    contour_levels = np.arange(lower_value, higher_value, interval)
    cmap = matplotlib.cm.get_cmap(cmap)
    return contour_levels, cmap


def colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, interval):

    colour_bar.set_ticks(np.arange(lower_tick, upper_tick, interval))
    colour_bar.set_ticklabels(np.arange(lower_tick, upper_tick, interval))
    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0, labelsize=8)

    return colour_bar

def plot_map(cubelist, variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval, subplot_columns, subplot_rows):

    crs_latlon = ccrs.PlateCarree()

    fig = plt.figure(figsize=(8, 8))

    cube_number = 0

    for i in cubelist:
        long_name = i.long_name
        ax = plt.subplot(subplot_columns, subplot_rows, cube_number+1, projection=crs_latlon)
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

    """Add colour bar."""
    if len(cubelist) > 20:
        colourbar_axis = fig.add_axes([0.20, 0.07, 0.60, 0.02])
    if len(cubelist) <=20:
        colourbar_axis = fig.add_axes([0.24, 0.42, 0.55, 0.015])

    colour_bar = plt.colorbar(contour_plot, colourbar_axis, orientation='horizontal')

    """Adjust ticks."""
    colour_bar = colour_bar_adjust_ticks(fig, contour_plot, colour_bar, lower_tick, upper_tick, tick_interval)

    if variable == 'pr':
        label = 'Precipitation (mm $\mathregular{day^{-1}}$)'
    if variable == 'hfls':
        label = 'Surface Upward Latent Heat Flux (W $\mathregular{m^{-2}}$)'
    if variable == 'evapotranspiration':
        label = 'Evaporation (mm $\mathregular{day^{-1}}$)'
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
    if variable == 'mrso_':
        label = 'Total Soil Moisture Content (mm)'
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
        label = 'Total Runoff Flux (mm $\mathregular{day^{-1}}$)'
    if variable == 'treeFrac':
        label = 'Tree Cover Fraction (%)'

    if variable == 'treeFrac':
        colour_bar.set_label(label, fontsize=8)
    if variable != 'treeFrac':
        colour_bar.set_label(input_time+" "+label, fontsize=9)

    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.47, top=0.95, wspace=0.45, hspace=0.5)
    print variable
    """Save the figure, close the plot and print an end statement."""
    if input_time in ['DJF', 'MAM', 'JJA', 'SON']:
        print season_index
        print "saving final figure"
        fig.savefig("Gridded_"+variable+"_"+str(season_index)+".png", bbox_inches='tight')
    if input_time in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:
        print month_index
        print "saving final figure"
        fig.savefig("Gridded_"+variable+"_"+str(month_index)+".png", bbox_inches='tight')
    plt.close()
    print "plot done"


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    list_of_reanalysis = ["gleam", "merra2"]
    #list_of_times = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'DJF', 'MAM', 'JJA', 'SON']
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CMCC-CM", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    model_type = "amip"


    list_of_times = ['MAM']
    #list_of_times = ['DJF', 'MAM', 'JJA', 'SON']

    variable = "hfls"
    lower_year = 1979
    upper_year = 2008
    lower_value = 80
    higher_value = 155
    value_interval = 5
    lower_tick = 80
    upper_tick = 160
    tick_interval = 10
    ensemble = "yes"
    cmap = "YlGnBu"
    subplot_columns = 4
    subplot_rows = 5

    # ------------------------------------------------------------------------------------------------------------------------------------------

    time_range_year = iris.Constraint(time=lambda cell: lower_year <= cell.point.year <= upper_year)

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # Loop through list of times

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

        """MODELS."""

        if len(list_of_models) == 0:
            pass
        else:

            """Extract the model file paths."""
            start_time = time.time()
            print list_of_models
            print model_type
            print variable
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

            """Slice each cube by year, month and spatial dimension."""
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

    # ------------------------------------------------------------------------------------------------------------------------------------------
            """ENSEMBLE MEAN."""

            if ensemble == "yes":

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

                """Compute the ensemble mean cube."""
                ensemble_mean_cube = sum(cubelist_ensemble) / float(len(cubelist_ensemble))
                ensemble_mean_cube.long_name = "Ensemble"

                """Append ensemble mean to the cubelist."""
                cubelist = np.append(cubelist, ensemble_mean_cube)

        # ------------------------------------------------------------------------------------------------------------------------------------------

        """REANALYSES."""

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
            for i in np.arange(0, len(cubelist_reanalysis)):
                p = multiprocessing.Process(target=slicing_reanalysis, args=(i, array))
                jobs.append(p)
                p.start()
            for process in jobs:
                process.join()
            cubelist_reanalysis = array.values()

            """Append reanalysis cubes to the cubelist."""
            cubelist = np.append(cubelist, cubelist_reanalysis)

            print cubelist

    #------------------------------------------------------------------------------------------------------------------------------------------

        """PLOT MAP."""

        start_time = time.time()
        """Define the contour levels and colour map for the plot."""
        contour_levels = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[0]
        cmap = contour_lev_colour_map(lower_value, higher_value, value_interval, cmap)[1]

        """Produce the plot."""
        plot_map(cubelist, variable, input_time, contour_levels, cmap, lower_tick, upper_tick, tick_interval, subplot_columns, subplot_rows)
        print time.time() - start_time, "seconds"
