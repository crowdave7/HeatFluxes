"""Import necessary modules for this code."""
import iris
import iris.analysis
import iris.coord_categorisation
import numpy as np
import os


def reanalysis_bowen(reanalysis_type, season_name):
    """Take the input variables, and find the paths to the relevant regridded sensible and latent heat flux model files."""
    """Plot these up as a spatial map for the given season."""
    """Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles"

    """Find the paths to the files containing the model data"""
    model_file_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            for j in reanalysis_type:
                if j in path and ('hfss' in path or 'hfls' in path):
                    model_file_paths = np.append(model_file_paths, path)

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    print model_file_paths

    """Load the sensible and latent heat flux data from each model into a cubelist."""

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    """Load the two cubes."""
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
            iris.coord_categorisation.add_season(i, 'time', name='clim_season')
            iris.coord_categorisation.add_season_year(i, 'time', name='season_year')

            """Aggregate the data by season and season year."""
            seasonal_means = i.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)

            """Constrain the data by the input season. Decapitalise the input variable."""
            constraint = iris.Constraint(clim_season=season_name.lower())
            reanalysis_data_season = seasonal_means.extract(constraint)

            """Take the mean over the cube."""
            reanalysis_data = reanalysis_data_season.collapsed('time', iris.analysis.MEAN)

            """Replace the cube in the cubelist with the new cube."""
            cubes[cube_id] = reanalysis_data

            """Add 1 to the count variable."""
            cube_id +=1

    """Select the numerator (latent heat flux)."""

    variable_name_0 = str(cubes[0].long_name)
    variable_name_1 = str(cubes[1].long_name)

    if "Latent" in variable_name_0 or "latent" in variable_name_0:
        numerator = cubes[0]
    elif "Latent" in variable_name_1 or "latent" in variable_name_1:
        numerator = cubes[1]

    """Calculate the denominator (latent + sensible heat flux)."""
    denominator = iris.analysis.maths.add(cubes[0], cubes[1])

    """Calculate the evaporative fraction (latent/(latent+sensible))."""
    evap_fraction = iris.analysis.maths.divide(numerator, denominator)

    if reanalysis_type == ["cfsr"]:
        evap_fraction.long_name = "CFSR"

    if reanalysis_type == ["doe"]:
        evap_fraction.long_name = "NCEP DOE-2"

    if reanalysis_type == ['erai']:
        evap_fraction.long_name = "ERA-Interim"

    if reanalysis_type == ['gleam']:
        evap_fraction.long_name = "GLEAM-LE"

    if reanalysis_type == ["jra"]:
        evap_fraction.long_name = "JRA-55"

    if reanalysis_type == ['merra2']:
        evap_fraction.long_name = "MERRA-2"

    if reanalysis_type == ["ncep"]:
        evap_fraction.long_name = "NCEP/NCAR"

    return evap_fraction

#reanalysis_bowen(["erai"], "SON")
