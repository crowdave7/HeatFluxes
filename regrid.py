#!/usr/bin/env python

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
    if variable == 'nrad':
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1NetRadiationModelFiles"
    if variable == 'vpd' or variable == 'hurs' or variable == 'tas' or variable == 'ws' or variable == 'swc_anom':
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles"

    if variable not in ["nrad", "evap_fraction", "vpd", "hurs", "tas", "ws", "swc_anom"]:
        #root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
        root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles"

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

    if variable in ['nrad', 'vpd', 'hurs', 'tas', 'ws', 'swc_anom']:
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

    if variable not in ['nrad', 'evap_fraction', 'vpd', 'hurs', 'tas', 'ws', 'swc_anom']:

        if root_directory == "/ouce-home/students/kebl4396/Paper1/Paper1CorrectedModelFiles":
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

        else:
            """Find the paths to the directories containing the model data"""
            print "hi1"
            directory_paths = []
            for root, directories, files in os.walk(root_directory):
                for i in directories:
                    path = os.path.join(root, i)
                    print path
                    print "hi2"
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

            model_file_paths_sorted = []

            for i in list_of_models:
                for j in model_file_paths:
                    if i in j:
                        model_file_paths_sorted = np.append(model_file_paths_sorted, j)

            model_file_paths = model_file_paths_sorted



    # """Find the paths to the directories containing the model data"""
    # directory_paths = []
    # print root_directory
    # for root, directories, files in os.walk(root_directory):
    #     for i in directories:
    #         path = os.path.join(root, i)
    #         for j in list_of_models:
    #             if j in path and model_type in path and ensemble in path:
    #                 directory_paths = np.append(directory_paths, path)
    #
    # """Find the model files and their absolute paths."""
    # model_file_paths = []
    # for i in directory_paths:
    #     files = os.listdir(i)
    #     for j in files:
    #         if variable in j:
    #             model_file_path = os.path.join(i, j)
    #             model_file_paths = np.append(model_file_paths, model_file_path)
    #
    # # model_file_paths = []
    # # for root, directories, files in os.walk(root_directory):
    # #     for i in files:
    # #         path = os.path.join(root, i)
    # #         for j in list_of_models:
    # #             if j in path and model_type in path and ensemble in path:
    # #                 model_file_paths = np.append(model_file_paths, path)

    model_file_paths = sorted(model_file_paths, key=lambda s: s.lower())

    print model_file_paths

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
    if variable == 'nrad':
        cubes = iris.load(model_file_paths)
    if variable == 'prveg':
        cubes = iris.load(model_file_paths)
    if variable == 'treeFrac':
        cubes = iris.load(model_file_paths)
    if variable == 'lai':
        cubes = iris.load(model_file_paths, 'leaf_area_index')
    if variable == 'mrros':
        cubes = iris.load(model_file_paths, 'surface_runoff_flux')
    if variable == 'mrro':
        cubes = iris.load(model_file_paths, 'runoff_flux')
        print "hi"
    if variable == 'mrso':
        cubes = iris.load(model_file_paths, 'soil_moisture_content')
    if variable == 'vpd':
        cubes = iris.load(model_file_paths)
    if variable == 'hurs':
        cubes = iris.load(model_file_paths)
    if variable == 'tas':
        cubes = iris.load(model_file_paths)
    if variable == 'ws':
        cubes = iris.load(model_file_paths)
    if variable == 'swc_anom':
        cubes = iris.load(model_file_paths)

    """Define the reanalysis data to regrid onto."""
    reanalysis_data = iris.load_cube("/ouce-home/data_not_backed_up/analysis/erai/1.0x1.0/daily/precip/nc/erai.totprecip.dt.1979-2013.nc")

    print cubes

    """For each cube,"""
    for cube in cubes:

        """Regrid the cube to the reanalysis data."""
        regridded_model_data = cube.regrid(reanalysis_data, iris.analysis.Linear())

        """Slice the regridded cube down to the African domain."""
        regridded_model_data = regridded_model_data.intersection(latitude=(-40, 40), longitude=(-30, 70))
        #
        # model_id = regridded_model_data.long_name

        """Extract the model ID for naming the netCDF file."""

        if variable in ['nrad', 'vpd', 'hurs', 'tas', 'ws', 'swc_anom']:
            model_id = regridded_model_data.long_name
        else:
            model_id = regridded_model_data.attributes['model_id']

        if model_id == 'ACCESS1.3':

            """Save the regridded cube as a netCDF file."""
            iris.save(regridded_model_data, "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles/"+variable+"_ACCESS1-3_"+model_type+"_regridded_Africa.nc")

        else:

            """Save the regridded cube as a netCDF file."""
            iris.save(regridded_model_data, "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles/"+variable+"_"+model_id+"_"+model_type+"_regridded.nc")

        """Print a statement to signify that regridding has finished for this model."""
        print model_id+" regridding done"

#regrid(["bcc-csm1-1-cd .m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"], "amip", "swc_anom")
#regrid(["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "evspsblveg")
#regrid(["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "tran")

#regrid(["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "pr")
#regrid(["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "hfls")
#regrid(["ACCESS1-0", "ACCESS1-3", "bcc-csm1-1_", "bcc-csm1-1-m", "BNU-ESM", "CanAM4", "CNRM-CM5", "CSIRO-Mk3-6-0", "GFDL-HIRAM-C180", "GFDL-HIRAM-C360", "GISS-E2-R", "HadGEM2-A", "inmcm4", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MPI-ESM-MR", "MRI-AGCM3-2H", "MRI-AGCM3-2S", "MRI-CGCM3", "NorESM1-M"], "amip", "nrad")
#regrid(["bcc-csm1-1_"], "amip", "nrad")
#regrid(["bcc-csm1-1/"], "amip", "lai")

#regrid(["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "vpd")

#regrid(["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "ws")
regrid(["CNRM-CM5/"], "amip", "mrsos")

#regrid(["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5", "CSIRO-Mk3-6-0/", "EC-EARTH", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR", "IPSL-CM5A-MR/", "IPSL-CM5B-LR", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"], "amip", "mrro")
