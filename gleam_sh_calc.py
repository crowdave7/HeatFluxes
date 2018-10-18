#!/usr/bin/env python

import iris
import iris.analysis
from iris.util import unify_time_units
import numpy as np
import netCDF4
from netCDF4 import num2date, date2num
import h5py
import iris.coords as icoords
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import os

print "hi"

#------------------------------------------------------------

# "CANOPY STORAGE CAPACITY EXTRACTION"
#
# list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]
#
# root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles"
#
# for m in list_of_models:
#
#     for char in '/':
#         model = m.replace(char,'')
#
#     model_name = str(model)
#
#     lai_file_path = []
#     for root, directories, files in os.walk(root_directory):
#         for i in files:
#             path = os.path.join(root, i)
#             #print path
#             if m == "bcc-csm1-1/":
#                 m = "bcc-csm1-1_"
#             for char in '/':
#                 m = m.replace(char,'')
#                 print m
#                 if m in path and 'amip' in path and 'lai' in path:
#                     lai_file_path = np.append(lai_file_path, path)
#
#     lai_cube = iris.load_cube(lai_file_path)
#
#     lai_cube = iris.analysis.maths.multiply(lai_cube, 0.2)
#
#     canopy_capacity_cube = iris.analysis.maths.add(lai_cube, 0.05)
#
#     canopy_capacity_cube.units = "mm day-1"
#
#     canopy_capacity_cube.long_name = model_name
#
#     iris.save(canopy_capacity_cube, "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles/cancap_"+model_name+"_amip_regridded.nc")

#-------------------------------------------------------------

# "SOIL MOISTURE ACCUMULATION / DEPLETION EXTRACTION FROM REGRIDDED MODEL FILES."
#
# list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
#
# root_directory = "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles"
#
# for m in list_of_models:
#
#     for char in '/':
#         model = m.replace(char,'')
#
#     model_name = str(model)
#
#     pr_file_path = []
#     for root, directories, files in os.walk(root_directory):
#         for i in files:
#             path = os.path.join(root, i)
#             if m == "bcc-csm1-1/":
#                 m = "bcc-csm1-1_"
#             for char in '/':
#                 m = m.replace(char,'')
#                 if m in path and 'amip' in path and 'pr_' in path:
#                     pr_file_path = np.append(pr_file_path, path)
#
#     hfls_file_path = []
#     for root, directories, files in os.walk(root_directory):
#         for i in files:
#             path = os.path.join(root, i)
#             if m == "bcc-csm1-1/":
#                 m = "bcc-csm1-1_"
#             for char in '/':
#                 m = m.replace(char,'')
#                 if m in path and 'amip' in path and 'hfls_' in path:
#                     hfls_file_path = np.append(hfls_file_path, path)
#
#     mrro_file_path = []
#     for root, directories, files in os.walk(root_directory):
#         for i in files:
#             path = os.path.join(root, i)
#             if m == "bcc-csm1-1/":
#                 m = "bcc-csm1-1_"
#             for char in '/':
#                 m = m.replace(char,'')
#                 if m in path and 'amip' in path and 'mrro_' in path:
#                     mrro_file_path = np.append(mrro_file_path, path)
#
#     print pr_file_path
#     print hfls_file_path
#     print mrro_file_path
#
#
#     pr_cube = iris.load_cube(pr_file_path)
#     pr_cube = iris.analysis.maths.multiply(pr_cube, 86400)
#
#     pr_cube.units = "mm day-1"
#
#     hfls_cube = iris.load_cube(hfls_file_path)
#     evap_cube = iris.analysis.maths.divide(hfls_cube, 28)
#
#     evap_cube.units = "mm day-1"
#
#     mrro_cube = iris.load_cube(mrro_file_path)
#     mrro_cube = iris.analysis.maths.multiply(mrro_cube, 86400)
#
#     mrro_cube.units = "mm day-1"
#
#     first_cube = iris.analysis.maths.subtract(pr_cube, evap_cube)
#
#     soil_accumulation_cube = iris.analysis.maths.subtract(first_cube, mrro_cube)
#
#     soil_accumulation_cube.long_name = model_name
#     soil_accumulation_cube.units = "mm day-1"
#
#     print soil_accumulation_cube
# #
#     iris.save(soil_accumulation_cube, "/ouce-home/students/kebl4396/Paper1/Paper1RegriddedModelFiles/sa_"+model_name+"_amip_regridded.nc")

#------------------------------------------------------------

"SOIL MOISTURE ACCUMULATION / DEPLETION EXTRACTION: MODEL FILES."

list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]

for m in list_of_models:

    for char in '/':
        model = m.replace(char,'')

    model_name = str(model)

    """Import the data."""
    root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    ensemble = "r1i1p1"
    model_type = "amip"

    """Find the paths to the directories containing the model data"""
    directory_paths = []
    for root, directories, files in os.walk(root_directory):
        for k in directories:
            path = os.path.join(root, k)
            if m in path and model_type in path and ensemble in path:
                directory_paths = np.append(directory_paths, path)

    print directory_paths

    """Find the model files to pr and their absolute paths."""
    pr_file_path = []
    for i in directory_paths:
        files = os.listdir(i)
        for j in files:
            if "pr_" in j:
                model_file_path = os.path.join(i, j)
                pr_file_path = np.append(pr_file_path, model_file_path)

    pr_file_path = sorted(pr_file_path, key=lambda s: s.lower())
    print pr_file_path

    if model_name in ["GFDL-HIRAM-C180", "GFDL-HIRAM-C360"]:
        pr_cube = iris.load(pr_file_path)[1]
    else:
        pr_cube = iris.load_cube(pr_file_path)

    """Find the model files to hfls and their absolute paths."""
    hfls_file_path = []
    for i in directory_paths:
        files = os.listdir(i)
        for j in files:
            if "hfls_" in j:
                model_file_path = os.path.join(i, j)
                hfls_file_path = np.append(hfls_file_path, model_file_path)

    hfls_file_path = sorted(hfls_file_path, key=lambda s: s.lower())
    print hfls_file_path

    if model_name in ["GFDL-HIRAM-C180", "GFDL-HIRAM-C360"]:
        hfls_cube = iris.load(hfls_file_path)[1]
    else:
        hfls_cube = iris.load_cube(hfls_file_path)

    """Find the model files to mrro and their absolute paths."""
    mrro_file_path = []
    for i in directory_paths:
        files = os.listdir(i)
        for j in files:
            if "mrro_" in j:
                model_file_path = os.path.join(i, j)
                mrro_file_path = np.append(mrro_file_path, model_file_path)

    mrro_file_path = sorted(mrro_file_path, key=lambda s: s.lower())
    print mrro_file_path

    if model_name in ["GFDL-HIRAM-C180", "GFDL-HIRAM-C360"]:
        mrro_cube = iris.load(mrro_file_path)[1]
    else:
        mrro_cube = iris.load_cube(mrro_file_path)

    print pr_file_path
    print hfls_file_path
    print mrro_file_path


    #pr_cube = iris.load_cube(pr_file_path)
    pr_cube = iris.analysis.maths.multiply(pr_cube, 86400)
    pr_cube.long_name = model_name
    pr_cube.units = "mm day-1"

    #hfls_cube = iris.load_cube(hfls_file_path)
    evap_cube = iris.analysis.maths.divide(hfls_cube, 28)
    evap_cube.long_name = model_name
    evap_cube.units = "mm day-1"

    #mrro_cube = iris.load_cube(mrro_file_path)
    mrro_cube = iris.analysis.maths.multiply(mrro_cube, 86400)
    mrro_cube.long_name = model_name
    mrro_cube.units = "mm day-1"

    print pr_cube
    print evap_cube
    print mrro_cube

    # iris.save(pr_cube, "/ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_pr_Amon_"+model_name+"_amip_r1i1p1.nc")
    # iris.save(evap_cube, "/ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_evap_Amon_"+model_name+"_amip_r1i1p1.nc")
    # iris.save(mrro_cube, "/ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_mrro_Amon_"+model_name+"_amip_r1i1p1.nc")

    # print pr_cube
    # print evap_cube
    # print mrro_cube

    first_cube = iris.analysis.maths.subtract(pr_cube, evap_cube)

    first_cube.long_name = model_name
    first_cube.units = "mm day-1"

    iris.save(first_cube, "/ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_pr_minus_evap_"+model_name+"_amip_r1i1p1.nc")
    iris.save(mrro_cube, "/ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_mrro_Amon_"+model_name+"_amip_r1i1p1.nc")



    # pr_cube.coord('latitude').attributes = {}
    # evap_cube.coord('latitude').attributes = {}
    # mrro_cube.coord('latitude').attributes = {}

    #second_cube = iris.analysis.maths.add(evap_cube, mrro_cube)


    # if first_cube.coord('latitude').points != mrro_cube.coord('latitude').points:
    #
    #     print "yes"

#
#     # first_cube.coord('latitude').attributes = {}
#

#
#     soil_accumulation_cube.long_name = model_name
#     soil_accumulation_cube.units = "mm day-1"
#
#     print soil_accumulation_cube
# #
    #iris.save(soil_accumulation_cube, "/ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/sa_"+model_name+"_amip.nc")

#------------------------------------------------------------

# "WIND SPEED MODEL EXTRACTION"
#
# """FUNCTIONS TO CALCULATE WIND SPEED USING PYTHAG."""
#
# def ws_data_func(u_data, v_data):
#     return np.sqrt(u_data**2 + v_data**2)
#
# def ws_units_func(u_cube, v_cube):
#     if u_cube.units != getattr(v_cube, 'units', u_cube.units):
#         raise ValueError("units do not match")
#     return u_cube.units
#
# #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
# #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
#
# list_of_models = ['GISS-E2-R/']
#
# for m in list_of_models:
#
#     #---------------------------------------------------------------------------------
#     "EXTRACT EASTWARD SURFACE WIND"
#
#     for char in '/':
#         model = m.replace(char,'')
#
#     model_name = str(model)
#
#     """Import the data."""
#     root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
#     ensemble = "r1i1p1"
#     model_type = "amip"
#
#     """Find the paths to the directories containing the model data"""
#     directory_paths = []
#     for root, directories, files in os.walk(root_directory):
#         for k in directories:
#             path = os.path.join(root, k)
#             if m in path and model_type in path and ensemble in path:
#                 directory_paths = np.append(directory_paths, path)
#
#     print directory_paths
#
#     """Find the model files to ua and their absolute paths."""
#     ua_file_path = []
#     for i in directory_paths:
#         files = os.listdir(i)
#         for j in files:
#             if "ua_" in j:
#                 model_file_path = os.path.join(i, j)
#                 ua_file_path = np.append(ua_file_path, model_file_path)
#
#     ua_file_path = sorted(ua_file_path, key=lambda s: s.lower())
#     print ua_file_path
#
#
#     if model_name in ["GFDL-HIRAM-C180", "GFDL-HIRAM-C360"]:
#         ua_cube = iris.load(ua_file_path)[1]
#     else:
#         ua_cube = iris.load_cube(ua_file_path)
#
#     uas_cube = ua_cube.extract(iris.Constraint(air_pressure=92500))
#
#     iris.save(uas_cube, "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles/atlas_uas_Amon_"+model_name+"_amip_r1i1p1.nc")

    # uas_cube = uas_cube.collapsed('time', iris.analysis.MEAN)
    #
    # qplt.contourf(uas_cube, 25)
    #
    # plt.gca().coastlines()
    #
    # plt.show()

    #---------------------------------------------------------------------------------
    # "EXTRACT NORTHWARD SURFACE WIND"
    #
    # for char in '/':
    #     model = m.replace(char,'')
    #
    # model_name = str(model)
    #
    # """Import the data."""
    # root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
    # ensemble = "r1i1p1"
    # model_type = "amip"
    #
    # """Find the paths to the directories containing the model data"""
    # directory_paths = []
    # for root, directories, files in os.walk(root_directory):
    #     for k in directories:
    #         path = os.path.join(root, k)
    #         if m in path and model_type in path and ensemble in path:
    #             directory_paths = np.append(directory_paths, path)
    #
    # print directory_paths
    #
    # """Find the model files to vas and their absolute paths."""
    # va_file_path = []
    # for i in directory_paths:
    #     files = os.listdir(i)
    #     for j in files:
    #         if "va_" in j:
    #             model_file_path = os.path.join(i, j)
    #             va_file_path = np.append(va_file_path, model_file_path)
    #
    # va_file_path = sorted(va_file_path, key=lambda s: s.lower())
    # print va_file_path
    #
    # if model_name in ["GFDL-HIRAM-C180", "GFDL-HIRAM-C360"]:
    #     va_cube = iris.load(va_file_path)[1]
    # else:
    #     va_cube = iris.load_cube(va_file_path)
    #
    # vas_cube = va_cube.extract(iris.Constraint(air_pressure=92500))

    # vas_cube = vas_cube.collapsed('time', iris.analysis.MEAN)
    #
    # qplt.contourf(vas_cube, 25)
    #
    # plt.gca().coastlines()
    #
    # plt.show()

    #iris.save(vas_cube, "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles/atlas_vas_Amon_"+model_name+"_amip_r1i1p1.nc")

#     #---------------------------------------------------------------------------------
#     """CALCULATE AND SAVE WIND SPEED."""
#
#     ws_ifunc = iris.analysis.maths.IFunc(ws_data_func, ws_units_func)
#     ws_cube = ws_ifunc(uas_cube, vas_cube, new_name='wind speed')
#
#     ws_cube.long_name = model_name
#     ws_cube.units = "m s-1"
#
#     print ws_cube
#
#     ws_cube = ws_cube.collapsed('time', iris.analysis.MEAN)
#
#     qplt.contourf(ws_cube, 25)
#
#     plt.gca().coastlines()
#
#     plt.show()
#
#     iris.save(ws_cube, "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles/atlas_ws_Amon_"+model_name+"_amip_r1i1p1.nc")
#
# #---------------------------------------------------------------------------------

# "SURFACE TEMP, RH AND VPD MODEL EXTRACTION"
#
#
# list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
#
# for m in list_of_models:
#
#     for char in '/':
#         model = m.replace(char,'')
#
#     model_name = str(model)
#
#     """Import the data."""
#     root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
#     ensemble = "r1i1p1"
#     model_type = "amip"
#
#     """Find the paths to the directories containing the model data"""
#     directory_paths = []
#     for root, directories, files in os.walk(root_directory):
#         for k in directories:
#             path = os.path.join(root, k)
#             if m in path and model_type in path and ensemble in path:
#                 directory_paths = np.append(directory_paths, path)
#
#     print directory_paths
#
#     """Find the model files to surface temp and their absolute paths."""
#     ta_file_path = []
#     for i in directory_paths:
#         files = os.listdir(i)
#         for j in files:
#             if "ta_" in j:
#                 model_file_path = os.path.join(i, j)
#                 ta_file_path = np.append(ta_file_path, model_file_path)
#
#     ta_file_path = sorted(ta_file_path, key=lambda s: s.lower())
#     print ta_file_path
#
#     if model_name in ["GFDL-HIRAM-C180", "GFDL-HIRAM-C360"]:
#         ta_cube = iris.load(ta_file_path)[1]
#     else:
#         ta_cube = iris.load_cube(ta_file_path)
#
#     tas_cube = ta_cube.extract(iris.Constraint(air_pressure=92500))
#
#     tas_cube.convert_units('celsius')
#
#     print tas_cube
#
#     tas_cube.long_name = model_name
#
#     print tas_cube
#
#     # iris.save(tas_cube, "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles/atlas_tas_Amon_"+model_name+"_amip_r1i1p1.nc")
#
#     # tas_cube = tas_cube.collapsed('time', iris.analysis.MEAN)
#     #
#     # qplt.contourf(tas_cube, 25)
#     #
#     # plt.gca().coastlines()
#     #
#     # plt.show()
#
#     numerator = iris.analysis.maths.multiply(tas_cube, 17.27)
#
#     denominator = iris.analysis.maths.add(tas_cube, 237.3)
#
#     second_half = iris.analysis.maths.divide(numerator, denominator)
#
#     exp = iris.analysis.maths.exp(second_half)
#
#     svp_cube = iris.analysis.maths.multiply(exp, 6.1078)
#
#     svp_cube.long_name = model_name
#
#     svp_cube.units = "hPa"
#
#     print svp_cube
#     #
#     # iris.save(svp_cube, "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles/atlas_svp_Amon_"+model_name+"_amip_r1i1p1.nc")
#
#     # svp_cube = svp_cube.collapsed('time', iris.analysis.MEAN)
#     #
#     # qplt.contourf(svp_cube, 25)
#     #
#     # plt.gca().coastlines()
#     #
#     # plt.show()
#
#     """Find the model files to RH and their absolute paths."""
#     hur_file_path = []
#     for i in directory_paths:
#         files = os.listdir(i)
#         for j in files:
#             if "hur_" in j:
#                 model_file_path = os.path.join(i, j)
#                 hur_file_path = np.append(hur_file_path, model_file_path)
#
#     hur_file_path = sorted(hur_file_path, key=lambda s: s.lower())
#     print hur_file_path
#
#     if model_name in ["GFDL-HIRAM-C180", "GFDL-HIRAM-C360"]:
#         hurs_cube = iris.load(hur_file_path)[1]
#         hurs_cube = hurs_cube.extract(iris.Constraint(air_pressure=92500))
#     else:
#         hurs_cube = iris.load_cube(hur_file_path, iris.Constraint(air_pressure=92500))
#
#     print hurs_cube
#
#     hurs_cube.long_name = model_name
#
#     print hurs_cube
#     #
#     # iris.save(hurs_cube, "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles/atlas_hurs_Amon_"+model_name+"_amip_r1i1p1.nc")
#
#     # hurs_cube = hurs_cube.collapsed('time', iris.analysis.MEAN)
#
#     # qplt.contourf(hurs_cube, 25)
#     #
#     # plt.gca().coastlines()
#     #
#     # plt.show()
#
#     # e/es x 100 = RH
#     #
#     # rh/100 = e/es
#     #
#     # rh/100 x es = e
#
#     rh_div_100 = iris.analysis.maths.divide(hurs_cube, 100)
#
#     avp_cube = iris.analysis.maths.multiply(rh_div_100, svp_cube)
#
#     avp_cube.long_name = "Actual Vapour Pressure"
#
#     avp_cube.units = "hPa"
#
#     print avp_cube
#
#     # iris.save(avp_cube, "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles/atlas_avp_Amon_"+model_name+"_amip_r1i1p1.nc")
#
#     # avp_cube = avp_cube.collapsed('time', iris.analysis.MEAN)
#     #
#     # qplt.contourf(avp_cube, 25)
#     #
#     # plt.gca().coastlines()
#     #
#     # plt.show()
#
#     vpd_cube = avp_cube - svp_cube
#
#     #vpd_cube = iris.analysis.maths.multiply(vpd_cube, -1)
#
#     vpd_cube = iris.analysis.maths.divide(vpd_cube, 10)
#
#     vpd_cube.long_name = model_name
#
#     avp_cube.units = "kPa"
#
#     print vpd_cube
#
#     iris.save(vpd_cube, "/ouce-home/students/kebl4396/Paper1/Paper1VPDModelFiles/atlas_vpd_Amon_"+model_name+"_amip_r1i1p1.nc")
#
#     # vpd_cube = vpd_cube.collapsed('time', iris.analysis.MEAN)
#     #
#     # qplt.contourf(vpd_cube, 25)
#     #
#     # plt.gca().coastlines()
#     #
#     # plt.show()

#------------------------------------------------------------

#"NET RAD MODEL EXTRACTION"

# list_of_models = ["BNU-ESM/", "CanAM4/"]
#
# for m in list_of_models:
#
#     for char in '/':
#         model = m.replace(char,'')
#
#     model_name = str(model)
#
#     """Import the data."""
#     root_directory = "/ouce-home/data_not_backed_up/model/cmip5"
#     ensemble = "r1i1p1"
#     model_type = "amip"
#
#     """Find the paths to the directories containing the model data"""
#     directory_paths = []
#     for root, directories, files in os.walk(root_directory):
#         for k in directories:
#             path = os.path.join(root, k)
#             if m in path and model_type in path and ensemble in path:
#                 directory_paths = np.append(directory_paths, path)
#
#     print directory_paths
#
#     """Find the model files and their absolute paths."""
#     rlds_file_path = []
#     for i in directory_paths:
#         files = os.listdir(i)
#         for j in files:
#             if "rlds" in j:
#                 model_file_path = os.path.join(i, j)
#                 rlds_file_path = np.append(rlds_file_path, model_file_path)
#
#     rlds_file_path = sorted(rlds_file_path, key=lambda s: s.lower())
#     print rlds_file_path
#
#
#     """Find the model files and their absolute paths."""
#     rsds_file_path = []
#     for i in directory_paths:
#         files = os.listdir(i)
#         for j in files:
#             if "rsds" in j:
#                 model_file_path = os.path.join(i, j)
#                 rsds_file_path = np.append(rsds_file_path, model_file_path)
#
#     rsds_file_path = sorted(rsds_file_path, key=lambda s: s.lower())
#     print rsds_file_path
#
#     """Find the model files and their absolute paths."""
#     rlus_file_path = []
#     for i in directory_paths:
#         files = os.listdir(i)
#         for j in files:
#             if "rlus" in j:
#                 model_file_path = os.path.join(i, j)
#                 rlus_file_path = np.append(rlus_file_path, model_file_path)
#
#     rlus_file_path = sorted(rlus_file_path, key=lambda s: s.lower())
#     print rlus_file_path
#
#
#     """Find the model files and their absolute paths."""
#     rsus_file_path = []
#     for i in directory_paths:
#         files = os.listdir(i)
#         for j in files:
#             if "rsus" in j:
#                 model_file_path = os.path.join(i, j)
#                 rsus_file_path = np.append(rsus_file_path, model_file_path)
#
#     rsus_file_path = sorted(rsus_file_path, key=lambda s: s.lower())
#     print rsus_file_path
#
#     print "hi3"
#     print str(m)
#
#     lw_cube_down = iris.load_cube(rlds_file_path, "surface_downwelling_longwave_flux_in_air")
#     sw_cube_down = iris.load_cube(rsds_file_path, "surface_downwelling_shortwave_flux_in_air")
#     lw_cube_up = iris.load_cube(rlus_file_path, "surface_upwelling_longwave_flux_in_air")
#     sw_cube_up = iris.load_cube(rsus_file_path, "surface_upwelling_shortwave_flux_in_air")
#
#     downwelling_cube = iris.analysis.maths.add(lw_cube_down, sw_cube_down)
#
#     print downwelling_cube.shape
#
#     upwelling_cube = iris.analysis.maths.add(lw_cube_up, sw_cube_up)
#
#     print upwelling_cube.shape
#
#     net_downward_rad_cube = iris.analysis.maths.subtract(downwelling_cube, upwelling_cube)
#
#     print net_downward_rad_cube.shape
#
#     net_downward_rad_cube.long_name = model_name
#
#     iris.save(net_downward_rad_cube, "/ouce-home/students/kebl4396/atlas_nrad_Amon_"+model_name+"_amip_r1i1p1_197901-200812.nc")
#
#     print net_downward_rad_cube

    # if str(m) == 'bcc-csm1-1/':
    #     model_name = "bcc-csm1-1"
    #

    #
    # print net_downward_rad_cube
    #
    # iris.save(net_downward_rad_cube, "/ouce-home/students/kebl4396/Paper1/Paper1NetRadiationModelFiles/atlas_nrad_Amon_"+model_name+"_amip_r1i1p1_197901-200812.nc")

    #"/ouce-home/students/kebl4396/Paper1/atlas_nrad_Amon_"+str(i)+"_amip_r1i1p1_197901-200812.nc"


    # net_rad_cube = net_rad_cube.collapsed('time', iris.analysis.MEAN)
    #
    # qplt.contourf(net_rad_cube, 25)
    #
    # # Add coastlines to the map created by contourf.
    # plt.gca().coastlines()
    #
    # plt.show()


#------------------------------------------------------------

#"GLEAM NET RAD EXTRACTION"

# SAVE SW DOWNWELLING

# for i in list_of_yearmonths:
#     cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc", "surface_downwelling_shortwave_flux")
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/cube/downwelling/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc")
#     print "done "+str(i)+""
# #

# SAVE LW DOWNWELLING

# for i in list_of_yearmonths:
#     cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/srb_rel3.1_longwave_monthly_"+str(i)+".nc", "surface_downwelling_longwave_flux")
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/cube/downwelling/srb_rel3.1_longwave_monthly_"+str(i)+".nc")
#     print "done "+str(i)+""

#SAVE SW UPWELLING

# for i in list_of_yearmonths:
#     cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc", "surface_upwelling_shortwave_flux")
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/cube/upwelling/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc")
#     print "done "+str(i)+""
#
# #SAVE LW UPWELLING
#
# for i in list_of_yearmonths:
#     cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/srb_rel3.1_longwave_monthly_"+str(i)+".nc", "surface_upwelling_longwave_flux")
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/cube/upwelling/srb_rel3.1_longwave_monthly_"+str(i)+".nc")
#     print "done "+str(i)+""

# cube1 = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/cube/upwelling/srb_rel3.1_longwave_monthly.nc")
# cube2 = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/cube/upwelling/srb_rel3.0_shortwave_monthly_utc.nc")
# print cube1
# print cube2

#LOAD LWR DOWNWELLING AND SWR DOWNWELLING

# lw_cube_down = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/srb_rel3.1_longwave_monthly_downwelling.nc")
# #print lw_cube_down
# sw_cube_down = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/srb_rel3.0_shortwave_monthly_utc_downwelling.nc")
# #print sw_cube_down
#
# # LOAD LWR UPWELLING AND SWR UPWELLING
#
# lw_cube_up = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/srb_rel3.1_longwave_monthly_upwelling.nc")
# #print lw_cube_up
# sw_cube_up = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/srb_rel3.0_shortwave_monthly_utc_upwelling.nc")
# #print sw_cube_up
#
# # ADD LWRDOWN + SWRDOWN TO GET TOTAL DOWNWELLING
#
# downwelling_cube = iris.analysis.maths.add(lw_cube_down, sw_cube_down)
#
#
# # ADD LWRUP + SWRUP TO GET TOTAL UPWELLING
# upwelling_cube = iris.analysis.maths.add(lw_cube_up, sw_cube_up)
#
# # DOWNWELLING MINUS UPWELLING GIVES NET DOWNWARD SOLAR Radiation
#
# net_downward_rad_cube = iris.analysis.maths.subtract(downwelling_cube, upwelling_cube)
# net_downward_rad_cube.long_name = "Net Downward Solar Radiation"
# print net_downward_rad_cube
# iris.save(net_downward_rad_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/srb_net_downward_rad.nc")




# CONSTRAIN YEARS IN HFLS CUBE
# #
# cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/evap/E_GLEAM_v3.2a_monthly.nc")
# time_range = iris.Constraint(time=lambda cell: 1980 <= cell.point.year <= 2008)
# with iris.FUTURE.context(cell_datetime_objects=True):
#     cube = cube.extract(time_range)
#     cube = iris.analysis.maths.multiply(cube, 28)
#     cube.long_name = "Surface Upward Latent Heat Flux"
#     cube.units = "W m-2"
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/evap/LH_GLEAM_v3.2a_monthly.nc")

# cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/tran_merra2.nc")
#
# with iris.FUTURE.context(cell_datetime_objects=True):
#     cube = iris.analysis.maths.divide(cube, 28)
#
#     cube.long_name = "Transpiration"
#
#     cube.units = "mm day-1"
#
#     print cube
#
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/tran_merra2_1.nc")

# # #LOAD NET RAD CUBE
#
# srb_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/nrad_gewex")
# print srb_cube
# #
#
# #LOAD HFLS CUBE
# #
# gleam_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/evap/E_GLEAM_v3.2a_monthly_1984_2007.nc")
# print gleam_cube
#
# #REGRID HFLS CUBE ONTO NET RAD CUBE
#
# srb_cube.coord('longitude').guess_bounds()
# srb_cube.coord('latitude').guess_bounds()
#
# gleam_cube.coord('longitude').guess_bounds()
# gleam_cube.coord('latitude').guess_bounds()
#
# #
# regridded_gleam_cube = gleam_cube.regrid(srb_cube, iris.analysis.AreaWeighted())
# print regridded_gleam_cube
# iris.save(regridded_gleam_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/evap/E_GLEAM_v3.2a_monthly_1984_2007_regridded.nc")


# # LOAD REGRIDDED HFLS CUBE AND NET RAD CUBE
#
# gleam_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/LH_GLEAM_v3.2a_monthly_1984_2007_regridded.nc")
# print gleam_cube.coord('time')
# srb_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/srb_net_downward_rad.nc")
# print srb_cube.coord('time')
#
# #REMOVE TIME COORD IN HFLS CUBE AND REPLACE WITH COORD IN NET RAD CUBE
#
# srb_coord = srb_cube.coord('time')
# gleam_cube.remove_coord('time')
# gleam_cube.add_dim_coord(srb_coord, 0)
#
# # SAVE NEW HFLS CUBE
#
# iris.save(gleam_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/LH_GLEAM_v3.2a_monthly_1984_2007_regridded_correct_coords.nc")
#
# # RELOAD HFLS AND SRB cubes

# srb_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/srb_net_downward_rad.nc")
#
# gleam_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/LH_GLEAM_v3.2a_monthly_1984_2007_regridded_correct_coords.nc")
# print gleam_cube.coord('time')
# print srb_cube.coord('time')
#
# #HFSS = NET DOWNWARD RAD - HFLS
#
# sensible_cube = iris.analysis.maths.subtract(srb_cube, gleam_cube)
#
# sensible_cube.long_name = "Surface Upward Sensible Heat Flux"
# sensible_cube.units = "W m-2"
# #
# # SAVE HFSS CUBE (180x360 GRID)
#
# iris.save(sensible_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/SH_GLEAM_estimate_v3.2a_monthly_1984_2007.nc")

# LOAD SENSIBLE CUBE
# LOAD LATENT CUBE
#
# nt_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/erai_nt.nc")
# ns_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/erai_ns.nc")
#
# print nt_cube
# print ns_cube
#
# nt_cube_2 = iris.analysis.maths.divide(nt_cube, -43200)
# ns_cube_2 = iris.analysis.maths.divide(ns_cube, 43200)
#
# iris.save(nt_cube_2, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/nt_cube_2.nc")
#
# iris.save(ns_cube_2, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/ns_cube_2.nc")
#
# nrad_cube = iris.analysis.maths.add(nt_cube, ns_cube)
# #
# # print nrad_cube
# #
# iris.save(nrad_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/nrad_cube.nc")

# cfsr_downward_lw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/cfsr_downward_lw.nc")
#
# cfsr_downward_sw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/cfsr_downward_sw.nc")
#
# cfsr_upward_lw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/cfsr_upward_lw.nc")
#
# cfsr_upward_sw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/cfsr_upward_sw.nc")
#
# print cfsr_downward_lw[1].data
# print cfsr_downward_sw[0].data
# print cfsr_upward_lw[0].data
# print cfsr_upward_sw[1].data
#
# iris.save(cfsr_downward_lw[1], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/cfsr_downward_lw_1.nc")
# iris.save(cfsr_downward_sw[0], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/cfsr_downward_sw_1.nc")
# iris.save(cfsr_upward_lw[0], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/cfsr_upward_lw_1.nc")
# iris.save(cfsr_upward_sw[1], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/cfsr_upward_sw_1.nc")


# print nt_cube
# print ns_cube
# #
# nt_cube_2 = iris.analysis.maths.divide(nt_cube, 86400)
# ns_cube_2 = iris.analysis.maths.divide(ns_cube, 86400)
#
# print nt_cube_2
# print ns_cube_2
#
# iris.save(nt_cube_2, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/nt_cube_2.nc")
#
# iris.save(ns_cube_2, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/ns_cube_2.nc")

# nrad_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/nrad_cube.nc")
#
# nrad_cube_2 = iris.analysis.maths.divide(nrad_cube, 86400)
#
# iris.save(nrad_cube_2, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/ERAI_data/net_rad/nrad_cube_2.nc")
# print hfss_cube
# print hfls_cube

# merra_sw_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/MERRA_data/net_rad/merra2_net_shortwave.nc")
# merra_lw_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/MERRA_data/net_rad/merra2_net_longwave.nc")
#
# merra_nrad_cube = iris.analysis.maths.add(merra_sw_cube, merra_lw_cube)

# ncepdoe_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/nrad_ncepdoe.nc")
# print ncepdoe_cube
# print ncepdoe_cube.long_name
#
# merra_nrad_cube.long_name = "Net Downward Solar Radiation"
#
# merra_nrad_cube.units = "W m-2"
#
# iris.save(merra_nrad_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/MERRA_data/net_rad/nrad_merra2.nc")

# gewex_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/nrad_merra2.nc")
#
# print gewex_cube.long_name

# jra_downward_lw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/JRA_data/JRA_downward_lw.nc")
#
# jra_downward_sw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/JRA_data/JRA_downward_sw.nc")
#
# jra_upward_lw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/JRA_data/JRA_upward_lw.nc")
#
# jra_upward_sw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/JRA_data/JRA_upward_sw.nc")
#
# iris.save(jra_downward_lw[0], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/JRA_data/JRA_downward_lw_1.nc")
# iris.save(jra_downward_sw[0], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/JRA_data/JRA_downward_sw_1.nc")
# iris.save(jra_upward_lw[1], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/JRA_data/JRA_upward_lw_1.nc")
# iris.save(jra_upward_sw[1], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/JRA_data/JRA_upward_sw_1.nc")
#
# print jra_downward_lw
# print jra_downward_sw
# print jra_upward_lw
# print jra_upward_sw

# cfsr_downward_lw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/CFSR_downward_lw.nc")
#
# cfsr_downward_sw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/CFSR_downward_sw.nc")
#
# cfsr_upward_lw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/CFSR_upward_lw.nc")
#
# cfsr_upward_sw = iris.load("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/CFSR_upward_sw.nc")
#
# iris.save(cfsr_downward_lw[1], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/CFSR_downward_lw_1.nc")
# iris.save(cfsr_downward_sw[0], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/CFSR_downward_sw_1.nc")
# iris.save(cfsr_upward_lw[1], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/CFSR_upward_lw_1.nc")
# iris.save(cfsr_upward_sw[1], "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/CFSR_data/net_rad/CFSR_upward_sw_1.nc")

# print cfsr_downward_lw
# print cfsr_downward_sw
# print cfsr_upward_lw
# print cfsr_upward_sw

#gleam_cube.coord('time').points = gleam_times

#print gleam_cube.coord('time')

#gleam_cube.remove_coord('time')

#gleam_cube.add_dim_coord(icoords.DimCoord(gleam_times), 0)

# era5_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1Code/pr_era5.nc")
#
# era5_cube_2 = iris.analysis.maths.multiply(era5_cube, 1000)
#
# iris.save(era5_cube_2, "/ouce-home/students/kebl4396/Paper1/Paper1Code/pr_era5_2.nc")
#
# print era5_cube_2


# list_of_yearmonths = []
# for i in np.arange(1984, 2008, 1):
#     for j in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
#         yearmonth = str(i) + j
#         list_of_yearmonths = np.append(list_of_yearmonths, yearmonth)
# print list_of_yearmonths

# for i in list_of_yearmonths:
#     cube_generator = iris.fileformats.netcdf.load_cubes("https://opendap.larc.nasa.gov/opendap/hyrax/SRB/GLW/SRB_REL3.1_LW_MONTHLY_NC/srb_rel3.1_longwave_monthly_"+str(i)+".nc")
#     cube = list(cube_generator)
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/srb_rel3.1_longwave_monthly_"+str(i)+".nc")
#     print "done "+str(i)+""
#
# list_of_yearmonths = []
# for i in np.arange(1984, 2008, 1):
#     for j in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
#         yearmonth = str(i) + j
#         list_of_yearmonths = np.append(list_of_yearmonths, yearmonth)
# print list_of_yearmonths

# for i in list_of_yearmonths:
#     cube_generator = iris.fileformats.netcdf.load_cubes("https://opendap.larc.nasa.gov/opendap/hyrax/SRB/GSW/SRB_REL3.0_SW_MONTHLY_UTC_NC/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc")
#     cube = list(cube_generator)
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc")
#     print "done "+str(i)+""

#------------------------------------------------------------

# """MODIS EVAPORATION EXTRACTION"""
#
# root_directory = "/ouce-home/data_not_backed_up/satellite/modis/8kmx8km/monthly/et/nc"
# file_paths = []
#
# for root, directories, files in os.walk(root_directory):
#     for i in files:
#         path = os.path.join(root, i)
#         file_paths = np.append(file_paths, path)
#
# #print file_paths
#
# file_paths_sorted = sorted(file_paths, key=lambda i: os.path.splitext(os.path.basename(i))[0])
#
# for i in file_paths_sorted:
#     cube = iris.load(i)
#     evap_cube = cube[2]
#     file_name = os.path.splitext(os.path.basename(i))[0]
#     print file_name
#     iris.save(evap_cube, "/ouce-home/projects/land_surface_climate/modis/"+file_name+".nc")
