"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.coord_categorisation
import numpy as np
import os
import iris.unit


def map_cube(list_of_models, model_type, variable, season_name):
    """Take the input variables, and find the paths to the relevant regridded model files."""
    """Produce the cube for the files."""
    """Central African domain."""

    """Import the data."""
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

    if variable == 'nrad':
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

    else:

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

    """If variable is evspsbl_, convert variable back to evspsbl"""
    if variable == 'evspsbl_':
        variable = 'evspsbl'

    """Load the data into a cubelist."""

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    cubes = []

    if variable == 'hfls':
        name = 'surface_upward_latent_heat_flux'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'hfss':
        name = 'surface_upward_sensible_heat_flux'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'pr':
        name = 'precipitation_flux'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'nrad':
        for i in model_file_paths:
            cube = iris.load_cube(i)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = i.long_name
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print len(times)
                print model_id
                print times[0]
                print times[-1]
                count +=1

    if variable == 'mrsos':
        name = 'moisture_content_of_soil_layer'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'tran':
        name = 'transpiration_flux'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'evspsblsoi':
        name = 'water_evaporation_flux_from_soil'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'evspsbl':
        name = 'water_evaporation_flux'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    if variable == 'evspsblveg':
        name = 'water_evaporation_flux_from_canopy'
        for i in model_file_paths:
            cube = iris.load_cube(i, name)
            cubes = np.append(cubes, cube)
        count = 0
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                variable_name = str(cubes[0].long_name)
                variable_units = str(cubes[0].units)
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

    """Take the mean over time for each cube in the cubelist."""
    cube_id = 0

    """For each cube,"""
    for regridded_model_data in cubes:

        """Select the model ID"""
        if variable == 'nrad':
            model_id = regridded_model_data.long_name
        else:
            model_id = cubes[cube_id].attributes['model_id']

        print model_id

        """If the variable is precipitation or transpiration, multiply model cubes by 86400"""
        if variable == 'pr' or variable == 'tran' or variable == 'evspsblsoi' or variable == 'evspsbl' or variable == 'evspsblveg':
            regridded_model_data = iris.analysis.maths.multiply(cubes[cube_id], 86400)

            """Reassign a model ID to the new multiplied cube."""
            regridded_model_data.long_name = model_id

        if variable != 'pr' or variable != 'tran' or variable != 'evspsblsoi' or variable != 'evspsbl':
            regridded_model_data.rename(model_id)

        """ If the input month is defined as the whole year,"""
        if season_name == 'Climatology':

            """Take the mean over the cube."""
            model_data = regridded_model_data.collapsed('time', iris.analysis.MEAN)

            """Replace the cube in the cubelist with the new cube."""
            cubes[cube_id] = model_data

            """Add 1 to the count variable."""
            cube_id +=1

        """ If the input month is defined as a season,"""
        if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

            """Create two coordinates on the cube to represent seasons and the season year."""
            iris.coord_categorisation.add_season(regridded_model_data, 'time', name='clim_season')
            iris.coord_categorisation.add_season_year(regridded_model_data, 'time', name='season_year')

            """Aggregate the data by season and season year."""
            seasonal_means = regridded_model_data.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

            """Constrain the data by the input season. Decapitalise the input variable."""
            constraint = iris.Constraint(clim_season=season_name.lower())
            regridded_model_data_season = seasonal_means.extract(constraint)

            """Take the mean over the cube."""
            model_data = regridded_model_data_season.collapsed('time', iris.analysis.MEAN)

            """Replace the cube in the cubelist with the new cube."""
            cubes[cube_id] = model_data

            """Add 1 to the count variable."""
            cube_id +=1

    print "Model cubes done"
    return cubes

#map_cube(["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "evspsblveg", "SON")

#map_cube(["ACCESS1-3"], "amip", "nrad", "SON")
