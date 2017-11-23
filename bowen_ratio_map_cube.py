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
import numpy as np
import os


def map_cube(list_of_models, model_type, season_name):
    """Take the input variables, and find the paths to the relevant regridded sensible and latent heat flux model files."""
    """Plot these up as a spatial map for the given season."""
    """Central African domain."""

    """Import the data."""
    root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

    """Find the paths to the files containing the model data"""
    model_file_paths = []
    for root, directories, files in os.walk(root_directory):
        for i in files:
            path = os.path.join(root, i)
            for j in list_of_models:
                if j in path and model_type in path and ensemble in path and ('hfss' in path or 'hfls' in path):
                    model_file_paths = np.append(model_file_paths, path)

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    print model_file_paths

    """Load the sensible and latent heat flux data from each model into a cubelist."""

    """Define a time range to constrain the years of the data."""
    time_range = iris.Constraint(time=lambda cell: 1979 <= cell.point.year <= 2008)

    """Load the data from the model file paths into a cube. Constrain the input years"""
    """Print the model ID, length of time dimension, and first and last model dates."""

    cubes_ef = []

    """ For each model compute the evaporative fraction in an iris cube and plot."""
    for j in list_of_models:

        """Find the relevant model paths, and load into a cubelist."""
        paths_for_this_model = [k for k in model_file_paths if j in k]
        cubes = iris.load(paths_for_this_model, ['surface_upward_latent_heat_flux', 'surface_upward_sensible_heat_flux'])
        count = 0
        print cubes

        """Constrain the years. Print the model ID, length of time dimension, and first and last model dates."""
        for i in cubes:
            with iris.FUTURE.context(cell_datetime_objects=True):
                cubes[count] = i.extract(time_range)
                time_points = cubes[count].coord('time').points
                times = cubes[count].coord('time').units.num2date(time_points)
                model_id = cubes[count].attributes['model_id']
                print model_id
                print len(times)
                print times[0]
                print times[-1]
                count +=1

        """Take the mean over time for each cube in the cubelist."""
        cube_id = 0

        """Print the model id from one of the cubes."""
        model_id = cubes[0].attributes['model_id']

        """For each cube (for the sensible and latent heat flux data for one model)"""
        for regridded_model_data in cubes:

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

        """Now compute the evaporative fraction for the one model."""

        """Select the numerator (latent heat flux)."""
        variable_name_0 = str(cubes[0].standard_name)
        variable_name_1 = str(cubes[1].standard_name)
        if "latent" in variable_name_0:
            numerator = cubes[0]
        elif "latent" in variable_name_1:
            numerator = cubes[1]

        """Calculate the denominator (latent + sensible heat flux)."""
        denominator = iris.analysis.maths.add(cubes[0], cubes[1])

        """Calculate the evaporative fraction (latent/(latent+sensible))."""
        evap_fraction = iris.analysis.maths.divide(numerator, denominator)

        cubes_ef = np.append(cubes_ef, evap_fraction)

    return cubes_ef

#map_cube(["ACCESS1-3"], "amip", "JJA")
