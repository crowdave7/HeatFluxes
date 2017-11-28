"""Import necessary modules for this code."""

import iris
import iris.analysis
import iris.coord_categorisation
import numpy as np
import os


def regrid(list_of_models, model_type, variable):
    """Take the input variables, and extract the model data from the server."""
    """Regrid the data to a common reanalysis grid, and save the data as netCDF."""

    """Import the data."""
    root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"

    """If variable is pr, distinguish between pr and precipitable water to find model files."""
    if variable == 'pr':
        variable = 'pr_'

    """If variable is evspsbl, distinguish between evspsbl and evspsblsoi to find model files."""
    if variable == 'evspsbl':
        variable = 'evspsbl_'

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

    """Load the data from the model file paths into a cube."""
    if variable == 'hfls':
        name = 'surface_upward_latent_heat_flux'
        cubes = iris.load(model_file_paths, name)
    if variable == 'hfss':
        name = 'surface_upward_sensible_heat_flux'
        cubes = iris.load(model_file_paths, name)
    if variable == 'pr':
        name = 'precipitation_flux'
        cubes = iris.load(model_file_paths, name)
    if variable == 'mrsos':
        name = 'moisture_content_of_soil_layer'
        cubes = iris.load(model_file_paths, name)
    if variable == 'tran':
        name = 'transpiration_flux'
        cubes = iris.load(model_file_paths, name)
    if variable == 'evspsblsoi':
        name = 'water_evaporation_flux_from_soil'
        cubes = iris.load(model_file_paths, name)
    if variable == 'evspsbl':
        name = 'water_evaporation_flux'
        cubes = iris.load(model_file_paths, name)
    if variable == 'evspsblveg':
        name = 'water_evaporation_flux_from_canopy'
        cubes = iris.load(model_file_paths, name)

    """Define the reanalysis data to regrid onto."""
    reanalysis_data = iris.load_cube("/ouce-home/data_not_backed_up/analysis/erai/1.0x1.0/daily/precip/nc/erai.totprecip.dt.1979-2013.nc")

    print cubes
    """For each cube,"""
    for cube in cubes:

        """Regrid the cube to the reanalysis data."""
        regridded_model_data = cube.regrid(reanalysis_data, iris.analysis.Linear())

        """Slice the regridded cube down to the African domain."""
        regridded_model_data = regridded_model_data.intersection(latitude=(-40, 40), longitude=(-30, 70))

        """Extract the model ID for naming the netCDF file."""
        model_id = regridded_model_data.attributes['model_id']
        if model_id == 'ACCESS1.3':

            """Save the regridded cube as a netCDF file."""
            iris.save(regridded_model_data, "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles/"+variable+"_ACCESS1-3_"+model_type+"_regridded_Africa.nc")

        else:

            """Save the regridded cube as a netCDF file."""
            iris.save(regridded_model_data, "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles/"+variable+"_"+model_id+"_"+model_type+"_regridded_Africa.nc")

        """Print a statement to signify that regridding has finished for this model."""
        print model_id+" regridding done"


#regrid(["ACCESS1-3", "bcc-csm1-1/", "BNU-ESM", "CanAM4", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C360", "GISS-E2-R/", "HadGEM2-A", "inmcm4", "MIROC5", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M/"], "amip", "evspsbl_")
regrid(["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "evspsblveg")
