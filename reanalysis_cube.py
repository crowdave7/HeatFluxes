"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.coord_categorisation
import numpy as np
import os
import iris.unit


def reanalysis(reanalysis_type, variable, season_name):
    """Take the input variables, and find the paths to the relevant regridded model files."""
    """Plot these up as a spatial map for the given season."""
    """Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles"

    """If variable is pr, distinguish between pr and precipitable water to find model files."""
    if variable == 'pr':
        variable = 'pr_'

    """Find the paths to the files containing the model data"""
    model_file_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            for j in reanalysis_type:
                if j in path and variable in path:
                    model_file_paths = np.append(model_file_paths, path)

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    print model_file_paths

    """If variable is pr_, convert variable back to pr"""
    if variable == 'pr_':
        variable = 'pr'

    """Load the data into a cubelist."""

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    """Load the cubes."""
    cubes = iris.load(model_file_paths)
    count = 0

    """For the cubes constrain the time range."""
    for i in cubes:
        with iris.FUTURE.context(cell_datetime_objects=True):
            cubes[count] = i.extract(time_range)
            time_points = cubes[count].coord('time').points
            times = cubes[count].coord('time').units.num2date(time_points)
            model_id = reanalysis_type
            print model_id
            print len(times)
            print times[0]
            print times[-1]
            count +=1

    """Take the mean over time for each cube in the cubelist."""
    cube_id = 0

    """For each cube (for the sensible and latent heat flux data for one reanalysis dataset)"""
    for i in cubes:

        """Slice the regridded cube down to the African domain."""
        cube = i.intersection(latitude=(-40, 40), longitude=(-30, 70))

        """Reminder of time points."""
        time_points = cube.coord('time').points
        times = cube.coord('time').units.num2date(time_points)

        """Plot up a map for the file."""

        """ If the input month is defined as the whole year,"""
        if season_name == 'Climatology':

            """Take the mean over the cube."""
            reanalysis_data = cube.collapsed('time', iris.analysis.MEAN)

            """Replace the cube in the cubelist with the new cube."""
            cubes[cube_id] = reanalysis_data

            """Add 1 to the count variable."""
            cube_id +=1

        """ If the input month is defined as a season,"""
        if season_name in ['DJF', 'MAM', 'JJA', 'SON']:

            """Create two coordinates on the cube to represent seasons and the season year."""
            iris.coord_categorisation.add_season(cube, 'time', name='clim_season')
            iris.coord_categorisation.add_season_year(cube, 'time', name='season_year')

            """Aggregate the data by season and season year."""
            seasonal_means = cube.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

            """Constrain the data by the input season. Decapitalise the input variable."""
            constraint = iris.Constraint(clim_season=season_name.lower())
            reanalysis_data_season = seasonal_means.extract(constraint)

            """Take the mean over the cube."""
            reanalysis_data = reanalysis_data_season.collapsed('time', iris.analysis.MEAN)

            """Replace the cube in the cubelist with the new cube."""
            cubes[cube_id] = reanalysis_data

            """Add 1 to the count variable."""
            cube_id +=1

    if reanalysis_type == ["cfsr"]:
        cubes[0].long_name = "CFSR"

    if reanalysis_type == ["doe"]:
        cubes[0].long_name = "NCEP DOE-2"

    if reanalysis_type == ['erai']:
        cubes[0].long_name = "ERA-Interim"

    if reanalysis_type == ['gleam']:
        cubes[0].long_name = "GLEAM-LE"

    if reanalysis_type == ["jra"]:
        cubes[0].long_name = "JRA-55"

    if reanalysis_type == ['merra2']:
        cubes[0].long_name = "MERRA-2"

    if reanalysis_type == ['mswep']:
        cubes[0].long_name = "GLEAM (MSWEP)"

    if reanalysis_type == ["ncep"]:
        cubes[0].long_name = "NCEP/NCAR"

    return cubes[0]

#reanalysis(["cfsr"], "pr", "SON")
