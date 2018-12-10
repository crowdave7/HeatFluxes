FIGURE 1

if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    #list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['pr']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'no'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    #list_prescribed_veg_models =


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 13.0
    y_tick_interval = 1.0
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'


FIGURE 2


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    #list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['evaporation']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'no'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    #list_prescribed_veg_models =


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 5.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'




FIGURE 3 gridded_map_factorised.py


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ['CanAM4/', 'CNRM-CM5/', 'bcc-csm1-1-m/', 'MRI-AGCM3-2S/', 'IPSL-CM5A-MR/', 'GISS-E2-R/', 'NorESM1-M/', 'GFDL-HIRAM-C360/', 'BNU-ESM/', 'MIROC5/', 'inmcm4/']

    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #'CNRM-CM5/', 'CanAM4/', 'bcc-csm1-1-m/', 'IPSL-CM5A-MR/', 'MRI-AGCM3-2S/', 'NorESM1-M/', 'MIROC5/', 'GFDL-HIRAM-C360/', 'GISS-E2-R/', 'BNU-ESM/', 'inmcm4/'
    #'CNRM-CM5/', 'CanAM4/', 'bcc-csm1-1-m/', 'GISS-E2-R/', 'IPSL-CM5A-MR/', 'MRI-AGCM3-2S/', 'NorESM1-M/', 'GFDL-HIRAM-C360/', 'BNU-ESM/', 'MIROC5/', 'inmcm4/'
    #'bcc-csm1-1-m/', 'GISS-E2-R/', 'CNRM-CM5/', 'GFDL-HIRAM-C360/', 'MRI-AGCM3-2S/', 'CanAM4/', 'NorESM1-M/', 'IPSL-CM5A-MR/', 'inmcm4/', 'BNU-ESM/', 'MIROC5/'
    #'CanAM4/', 'CNRM-CM5/', 'bcc-csm1-1-m/', 'MRI-AGCM3-2S/', 'IPSL-CM5A-MR/', 'GISS-E2-R/', 'NorESM1-M/', 'GFDL-HIRAM-C360/', 'BNU-ESM/', 'MIROC5/', 'inmcm4/'
    #list_of_models = []
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['jra', 'merra2']
    #list_of_times = ['DJF', 'MAM', 'JJA', 'SON', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    #list_of_times = ['Jun', 'Oct', 'Nov', 'Dec']
    #list_of_times = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    list_of_times = ['Mar', 'Jul', 'Nov', 'Jan']
    #list_of_times = ['DJF', 'MAM', 'JJA', 'SON']
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CMCC-CM", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    model_type = "amip"


    #list_of_times = ['DJF']
    #list_of_times = ['DJF', 'MAM', 'JJA', 'SON']

    variable = "evaporation"
    lower_year = 1979
    upper_year = 2008
    lower_value = 1.0
    higher_value = 5.5
    value_interval = 0.5
    lower_tick = 1.0
    upper_tick = 5.5
    tick_interval = 0.5
    ensemble = "yes"
    cmap = "YlGnBu"
    unit_plot = "mm day-1"
    subplot_columns = 4
    subplot_rows = 4
    composite_1 = 'yes'
    composite_2 = 'yes'
    list_of_models_composite_1 = ["CNRM-CM5/", "GISS-E2-R/"]
    list_of_models_composite_2 = ["BNU-ESM/", "MIROC5/"]






FIGURE 4

if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    #list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['pr']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    #list_prescribed_veg_models =


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 13.0
    y_tick_interval = 1.0
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'



FIGURE 5

if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    #list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['evaporation']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    #list_prescribed_veg_models =


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 5.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'


FIGURE 6 correlation_2.py


#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Import necessary modules for this code."""

import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.basemap import Basemap, maskoceans, interp
import numpy as np
from scipy.stats.stats import pearsonr
matplotlib.use('Agg')

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on', labelsize=7)
ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on', labelsize=7)

#Evaporation v LandFlux-EVAL rmse
#variable_y = [0.52513046, 1.0515319, 0.30871224, 0.24810952, 0.77596768, 0.46037268, 1.23658859, 0.50148155, 1.02871562, 0.71197554, 0.61166957, 0.62293586, 0.56419109, 0.07051655]

# Precip v Mean (CHIRPS/GPCC) RMSE
#variable_x = [1.05976145, 3.41556894, 1.45201831, 0.63344791, 0.82590173, 0.99572356, 0.87939634, 3.01854602, 2.81287744, 1.42418161, 2.92457451, 1.49827281, 0.31040471, 0.1935775]

#Change in Precipitation Nov - Mar
variable_x = [0.442, 0.497, 1.317, -0.366, 1.416, 0.609, -0.187, -0.618, 1.049, 0.811, 0.321, 0.433, 0.132, 0.697, 0.292]

# Change in Evaporation Nov-Mar
variable_y = [0.024, -0.012, -0.395, -0.575, -0.031, -0.488, -0.6, -0.314, -0.009, -0.381, -0.152, -0.287, -0.545, -0.025, -0.430]

list_of_models = ['bcc-csm1-1-m', 'BNU-ESM', 'CanAM4', 'CNRM-CM5', 'GFDL-HIRAM-C360', 'GISS-E2-R', 'inmcm4', 'IPSL-CM5A-MR', 'MIROC5', 'MRI-AGCM3-2S', 'NorESM1-M', 'GCM Ensemble', 'Best GCMs', 'Worst GCMs', 'MSWEP/GLEAM']

print len(variable_x)
print len(variable_y)
print len(list_of_models)

cmap_models = plt.get_cmap('rainbow')
colours_models = [cmap_models(i) for i in np.linspace(0, 1, 11)]

cmap_reanalysis = plt.get_cmap('brg')
colours_reanalysis = [cmap_reanalysis(i) for i in np.linspace(0, 1, len(variable_x))]

"""For each model,"""
for i in np.arange(0, len(variable_x), 1):

    """Select the line colour and add one to the line count variable."""

    if list_of_models[i] == "AGCM Ensemble":
        ax1.scatter(variable_x[i], variable_y[i], color='black', marker='o', s=16, label=list_of_models[i])

    elif list_of_models[i] == "Best AGCMs":
        ax1.scatter(variable_x[i], variable_y[i], color='forestgreen', marker='o', s=16, label=list_of_models[i])

    elif list_of_models[i] == "Worst AGCMs":
        ax1.scatter(variable_x[i], variable_y[i], color='saddlebrown', marker='o', s=16, label=list_of_models[i])

    elif list_of_models[i] == "MSWEP/GLEAM":
        ax1.scatter(variable_x[i], variable_y[i], color=colours_reanalysis[0], marker='o', s=16, label=list_of_models[i])

    else:

        """Add the scatter plot. Select dot colour, marker and label."""
        ax1.scatter(variable_x[i], variable_y[i], color=colours_models[i], marker='o', s=16, label=list_of_models[i])

    """Add a legend for each dot."""
    legend = plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.525), fontsize=7)

plt.ylabel(u'Δ Evaporation between March/November (mm $\mathregular{day^{-1}}$)', fontsize=7)
plt.ylim((-0.7, 0.1))
plt.yticks([-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1])

plt.xlabel(u'Δ Rainfall between March/November (mm $\mathregular{day^{-1}}$)', fontsize=7)
plt.xlim((-1.0, 2.0))
plt.xticks([-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0])


# x = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]
#
# y = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]
#
# ax1.plot(x, y, color = 'k', linestyle = '--', linewidth=1)

# Move left y-axis and bottim x-axis to centre, passing through (0,0)
#ax1.spines['left'].set_position('zero')
#ax1.spines['bottom'].set_position('zero')

# Eliminate upper and right axes
#ax1.spines['right'].set_color('none')
#ax1.spines['top'].set_color('none')

# Show ticks in the left and lower axes only
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

plt.axvline(0, color = 'black', linewidth=0.9)
plt.axhline(0, color = 'black', linewidth=0.9)

#plt.title("Congo Basin (14$^\circ$S - 4$^\circ$N, 18$^\circ$E - 29$^\circ$E)")

"""Save figure."""
fig.savefig("Change_precip_evap_Congo.png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)












FIGURE 7

if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['evspsblveg']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    #list_prescribed_veg_models =


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 3.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'


FIGURE 8


if seasonal_cycle_dv_array != []:
    line = ax1.plot(x_pos, seasonal_cycle_dv, zorder=4, linestyle='-', linewidth=3.0, color='firebrick', label = "Dynamic Vegetation")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("-")
    if include_legend == 'yes':
        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

if seasonal_cycle_pv_array != []:
    line = ax1.plot(x_pos, seasonal_cycle_pv, zorder=4, linestyle='-', linewidth=3.0, color='dodgerblue', label = "Prescribed Vegetation")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("-")
    if include_legend == 'yes':
        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)


plus



FIGURE 9


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['tran']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    #list_prescribed_veg_models =


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 3.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'


FIGURE 10

model_line_colours = [(0.5, 0.0, 1.0, 1.0), (0.096078431372549011, 0.80538091938883261, 0.8924005832479478, 1.0), (0.30000000000000004, 0.95105651629515353, 0.80901699437494745, 1.0), (0.50392156862745097, 0.99998102734872685, 0.70492554690614717, 1.0), (0.69999999999999996, 0.95105651629515364, 0.58778525229247314, 1.0), (0.90392156862745088, 0.80538091938883272, 0.45124405704532283, 1.0)]

plus


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['tran']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    #list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    list_prescribed_veg_models = []


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 3.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'


FIGURE 11


model_line_colours = [(0.30392156862745101, 0.30315267411304353, 0.98816547208125938, 1.0), (0.099999999999999978, 0.58778525229247314, 0.95105651629515353, 1.0), (1.0, 0.58778525229247325, 0.30901699437494745, 1.0), (1.0, 0.30315267411304364, 0.15339165487868545, 1.0), (1.0, 1.2246467991473532e-16, 6.123233995736766e-17, 1.0)]

plus


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    list_of_reanalysis = []
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['evspsblsoi']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    #list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    list_dynamic_veg_models = []
    #list_prescribed_veg_models = []


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 3.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'






FIGURE 12

Kahiu et al. plot


FIGURE 13


if seasonal_cycle_dv_array != []:
    line = ax1.plot(x_pos, seasonal_cycle_dv, zorder=4, linestyle='-', linewidth=3.0, color='firebrick', label = "Dynamic Vegetation")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("-")
    if include_legend == 'yes':
        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)

if seasonal_cycle_pv_array != []:
    line = ax1.plot(x_pos, seasonal_cycle_pv, zorder=4, linestyle='-', linewidth=3.0, color='dodgerblue', label = "Prescribed Vegetation")
    handles, labels = ax1.get_legend_handles_labels()
    handles[-1].set_linestyle("-")
    if include_legend == 'yes':
        legend = plt.legend(handles, labels, bbox_to_anchor=(1.03, 0.5), loc="center left", fontsize=9, handlelength=2.5)


plus


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    list_of_reanalysis = []
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['tran']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    #list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    #list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    list_dynamic_veg_models = []
    #list_prescribed_veg_models = []


    list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 3.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'no'
    ensemble = 'yes'
    reanalysis = 'no'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'




FIGURE 14


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['evspsblsoi']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'yes'

    list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    #list_prescribed_veg_models =


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 3.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'


FIGURE 15


    #flat bse
    model_line_colours = [(0.30392156862745101, 0.30315267411304353, 0.98816547208125938, 1.0), (0.099999999999999978, 0.58778525229247314, 0.95105651629515353, 1.0), (0.096078431372549011, 0.80538091938883261, 0.8924005832479478, 1.0), (0.69999999999999996, 0.95105651629515364, 0.58778525229247314, 1.0), (1.0, 0.58778525229247325, 0.30901699437494745, 1.0)]


plus


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    list_of_reanalysis = ['gleam']
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['evspsblsoi']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'no'

    #list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    list_dynamic_veg_models = []
    #list_prescribed_veg_models = []


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 3.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'


FIGURE 16


    #spikey bse
    #model_line_colours = [(0.5, 0.0, 1.0, 1.0), (0.30000000000000004, 0.95105651629515353, 0.80901699437494745, 1.0), (0.50392156862745097, 0.99998102734872685, 0.70492554690614717, 1.0), (0.90392156862745088, 0.80538091938883272, 0.45124405704532283, 1.0), (1.0, 0.30315267411304364, 0.15339165487868545, 1.0), (1.0, 1.2246467991473532e-16, 6.123233995736766e-17, 1.0)]


plus


if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------------------------------------------------
    # LIST OF INPUTS

    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CCSM4/", "CESM1-CAM5/", "CMCC-CM/", "CNRM-CM5/", "CSIRO-Mk3-6-0/", "EC-EARTH/", "FGOALS-g2/", "FGOALS-s2/", "GFDL-CM3/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/", "ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/",  "CSIRO-Mk3-6-0/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["ACCESS1-0/"]
    #list_of_models = []
    #evap from canopy
    #list_of_models = ["ACCESS1-0/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "HadGEM2-A/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #transpiration
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MPI-ESM-LR/", "MPI-ESM-MR/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #bare soil evaporation
    #list_of_models = ["ACCESS1-3/", "bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["CMCC-CM/"]
    #all for barchart
    #list_of_models = ["bcc-csm1-1-m/", "GFDL-HIRAM-C360/", "IPSL-CM5A-MR/", "MIROC5/", "NorESM1-M/"]

    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-LR/", "IPSL-CM5A-MR/", "IPSL-CM5B-LR/", "MIROC5/", "MRI-AGCM3-2H/", "MRI-AGCM3-2S/", "MRI-CGCM3/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]


    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["CNRM-CM5/", "GISS-E2-R/"]
    #list_of_models = ["bcc-csm1-1-m/"]
    #list_of_models = ["BNU-ESM/", "CNRM-CM5/"]
    #list_of_models = ["bcc-csm1-1-m/", "BNU-ESM/", "CanAM4/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/", "IPSL-CM5A-MR/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "MIROC5/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    list_of_models = ["bcc-csm1-1-m/",  "GFDL-HIRAM-C360/", "GISS-E2-R/", "IPSL-CM5A-MR/", "MRI-AGCM3-2S/", "NorESM1-M/"]
    #list_of_models = ["BNU-ESM/", "CanAM4/", "CNRM-CM5/","inmcm4/", "MIROC5/"]
    #list_of_models = []
    #list_of_models = []
    #list_of_reanalysis = ["era5", "jra", "merra2", "mswep"]
    #list_of_reanalysis = ["era5", "gleam", "jra", "merra2"]
    #list_of_reanalysis = ['era5', 'gleam', 'landfluxeval', 'modis']
    #list_of_reanalysis = ['chirps', 'era5', 'gpcc', 'mswep']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/", "CNRM-CM5/", "GFDL-HIRAM-C180/", "GFDL-HIRAM-C360/", "GISS-E2-R/", "inmcm4/"]
    #list_of_reanalysis = ['cfsr', 'era5', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_reanalysis = ["cfsr", "era5", "erai", "gleam", "jra", "merra2", "ncep-doe"]
    #list_of_reanalysis = ['gleam', 'merra2']
    #list_of_reanalysis = ['gewex', 'merra2']
    #list_of_reanalysis = ['cfsr', 'erai', 'gleam', 'jra', 'merra2', 'ncep-doe']
    #list_of_models = ["bcc-csm1-1/", "bcc-csm1-1-m/"]
    #list_of_reanalysis = ["cfsr", "era5", "erai", "jra", "merra2", "mswep", "ncep-doe"]
    #list_of_reanalysis = ['era5', 'gewex']
    #list_of_reanalysis = ['gleam', 'landfluxeval']
    #list_of_reanalysis = ['chirps', 'gpcc', 'mswep']
    list_of_reanalysis = []
    #list_of_reanalysis = ['chirps', 'mswep']

    model_type = "amip"

    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #partitioning one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation"]
    #normalisation one
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "evaporation", "pr"]
    #list_of_variables = ["hfls", "hfss", "nrad"]
    #list_of_variables = ["evspsblsoi", "tran", "evspsblveg", "hfls"]
    #list_of_variables = ['tran']
    #list_of_variables = ['evspsblsoi']
    #list_of_variables = ['pr', 'evspsblsoi', 'evspsblveg', 'tran', 'mrro']
    #list_of_variables = ['pr', 'evaporation', 'mrro']
    list_of_variables = ['evspsblsoi']


    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi"]
    #list_of_variables = ['evaporation']
    #list_of_variables = ['evaporation']
    #list_of_variables = ["evaporation", "evspsblveg", "tran", "evspsblsoi", "pr", "evaporation", "mrro", "sa", "nrad", "hfls", "hfss", 'tas', "hurs", "ws", 'tas', "vpd", "mrsos", "sa", "lai", "cancap"]
    #list_of_variables = ['tas', 'hurs', 'ws', 'vpd']
    subplot_multiple_variables = 'no'
    include_dynamic_prescribed = 'no'

    #list_dynamic_veg_models = ['CNRM-CM5/', 'GISS-E2-R/']
    list_prescribed_veg_models = ["BNU-ESM/", "MIROC5/"]
    list_dynamic_veg_models = []
    #list_prescribed_veg_models = []


    #list_dynamic_veg_models = ['bcc-csm1-1-m/', 'BNU-ESM/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    ##list_dynamic_veg_models = ['bcc-csm1-1-m/', 'GFDL-HIRAM-C360/', 'IPSL-CM5A-MR/', 'inmcm4/']
    #list_prescribed_veg_models = ['GISS-E2-R/', 'CNRM-CM5/']

    #list_dynamic_veg_models = ['BNU-ESM/']
    #list_prescribed_veg_models = ['MIROC5/', "MRI-AGCM3-2S/", "NorESM1-M/", 'CanAM4/', 'CNRM-CM5/']

    #list_prescribed_veg_models = []



    subplot_columns = 1
    subplot_rows = 2

    variables_to_add = []
    #variables_to_add = [['evaporation', 'mrro']]
    variables_to_subtract = []
    #variables_to_subtract = [['pr', 'evspsblsoi', 'evspsblveg', 'mrro']]
    #variables_to_subtract = [['pr', 'mrro']]
    #partitioning one
    #variables_to_subtract = [['evaporation', 'tran', 'evspsblveg', 'evspsblsoi']]
    #variables_to_subtract = [['pr', 'mrro']]
    divide_four_variables_by_fifth = 'no'

    # EXACT NUMBERS NOT PYTHON INDICES
    #fill_between_lines = [0, 3]
    #fill_between_lines = [3, 6]

    #fill_between_lines = [1, 4]
    fill_between_lines = []
    #variables_to_add = []

    #CONGO BASIN
    lower_lat = -14
    upper_lat = 4
    lower_lon = 18
    upper_lon = 29

    # CONGO RAINFOREST
    # lower_lat = -5
    # upper_lat = 4
    # lower_lon = 18
    # upper_lon = 29

    # CONGO BASIN SOUTH
    # lower_lat = -14
    # upper_lat = -5
    # lower_lon = 18
    # upper_lon = 29

    lower_year = 1979
    upper_year = 2008
    lower_y_lim = 0.0
    upper_y_lim = 3.5
    y_tick_interval = 0.5
    # lower_y_lim = -1
    # upper_y_lim = 5
    # y_tick_interval = 1.0
    depth = 0.85011389

    cmap_models = 'rainbow'
    cmap_reanalysis = 'brg'
    models = 'yes'
    ensemble = 'yes'
    reanalysis = 'yes'
    plot = 'yes'
    unit_plot = "mm day-1"
    bar_plot = 'no'
    bar_times = ['Mar']
    duplicate_seasonal_cycle = 'yes'
    bar_y_axis_title = "Evaporation (mm $\mathregular{day^{-1}}$)"
    #bar_y_axis_title = 'Normalised Evaporation (E/P)'
    #bar_y_axis_title = 'Heat Flux (W $\mathregular{m^{-2}}$)'
    #bar_y_axis_title = "Transpiration (mm $\mathregular{day^{-1}}$)"
    bar_colours = ['saddlebrown', 'dodgerblue', 'forestgreen']
    #bar_colours = ['lightseagreen', 'darkorange']
    bar_width = 0.5
    legend_in_plot = 'yes'
    include_legend = 'yes'
    rmse = 'no'
