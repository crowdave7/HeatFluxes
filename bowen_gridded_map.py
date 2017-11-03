"""Import necessary modules for this code."""
import bowen_ratio_ensemble_cube
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
import bowen_ratio_reanalysis_cube


def plot_bowen(list_of_models, model_type, season_name):
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

    """Set up an array to store evaporative fraction cubes."""
    evap_fraction_cubes = []

    """Define the contour levels for the input variables."""
    contour_levels = np.arange(0.4, 1.05, 0.05)

    """Define the colour map and the projection."""
    cmap = matplotlib.cm.get_cmap('YlGnBu')
    crs_latlon = ccrs.PlateCarree()

    """Plot the figure."""
    fig = plt.figure(figsize=(8, 8))
    #12, 6 for 5x5 subplot

    """ For each model compute the evaporative fraction in an iris cube and plot."""
    for j in list_of_models:

        """Find the relevant model paths, and load into a cubelist."""
        paths_for_this_model = [k for k in model_file_paths if j in k]
        cubes = iris.load(paths_for_this_model, ['surface_upward_latent_heat_flux', 'surface_upward_sensible_heat_flux'])
        count = 0

        """Constrain the years. Print the model ID, length of time dimension, and first and last model dates."""
        for i in cubes:
            print cubes
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

        evap_fraction.long_name = model_id

        evap_fraction_cubes = np.append(evap_fraction_cubes, evap_fraction)

    """Add the ensemble mean to the cubelist."""
    bowen_ensemble_cube_data = bowen_ratio_ensemble_cube.bowen_ensemble(list_of_models, model_type, season_name)
    evap_fraction_cubes = np.append(evap_fraction_cubes, bowen_ensemble_cube_data)

    """Add CFSR data to the cubelist."""
    cfsr_cube_data = bowen_ratio_reanalysis_cube.reanalysis_bowen(["cfsr"], season_name)
    evap_fraction_cubes = np.append(evap_fraction_cubes, cfsr_cube_data)

    """Add ERA-Interim data to the cubelist."""
    erai_cube_data = bowen_ratio_reanalysis_cube.reanalysis_bowen(["erai"], season_name)
    evap_fraction_cubes = np.append(evap_fraction_cubes, erai_cube_data)

    """Add GLEAM data to the cubelist."""
    gleam_cube_data = bowen_ratio_reanalysis_cube.reanalysis_bowen(["gleam"], season_name)
    evap_fraction_cubes = np.append(evap_fraction_cubes, gleam_cube_data)

    """Add JRA-55 data to the cubelist."""
    jra_cube_data = bowen_ratio_reanalysis_cube.reanalysis_bowen(["jra"], season_name)
    evap_fraction_cubes = np.append(evap_fraction_cubes, jra_cube_data)

    """Add MERRA-2 data to the cubelist."""
    merra2_cube_data = bowen_ratio_reanalysis_cube.reanalysis_bowen(["merra2"], season_name)
    evap_fraction_cubes = np.append(evap_fraction_cubes, merra2_cube_data)

    """Add NCEP DOE-2 data to the cubelist."""
    doe_cube_data = bowen_ratio_reanalysis_cube.reanalysis_bowen(["doe"], season_name)
    evap_fraction_cubes = np.append(evap_fraction_cubes, doe_cube_data)

    """Plot each model and the ensemble mean up."""
    model_number = 0

    for model_data in evap_fraction_cubes:

        """Plot the map using cartopy, and add map features."""
        ax = plt.subplot(6, 4, model_number+1, projection=crs_latlon)
        ax.set_extent([-22, 62, -24, 17], crs=crs_latlon)
        contour_plot = iplt.contourf(model_data, contour_levels, cmap=cmap, extend='both')

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

        """Add a title."""
        plt.title(model_data.long_name, fontsize=8)

        """Add 1 to the model number to loop through the next model."""
        model_number +=1
        print model_number

    colourbar_axis = fig.add_axes([0.20, 0.07, 0.60, 0.02])
    colour_bar = plt.colorbar(contour_plot, colourbar_axis, orientation='horizontal')

    colour_bar.set_ticks(np.arange(0.4, 1.05, 0.1))
    colour_bar.set_ticklabels(np.arange(0.4, 1.05, 0.1))

    colour_bar.ax.tick_params(axis=u'both', which=u'both', length=0)

    colour_bar.set_label(("Evaporative Fraction"), fontsize=10)

    fig.subplots_adjust(left=0.125, right=0.9, bottom=0.14, top=0.93, wspace=0.4, hspace=0.9)

    """Save the figure, close the plot and print an end statement."""
    print "saving final figure"
    fig.savefig("Evap_Frac_"+season_name+"_gridded.png", dpi=600, aspect='auto')
    plt.close()
    print " plot done"

#plot_bowen(["ACCESS1-3"], "amip", "SON")

plot_bowen(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CNRM-CM5/", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], "amip", "SON")
