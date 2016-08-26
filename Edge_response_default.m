%general main for edge response software

function Edge_response_default()

    %creates runvar structure variable that stores variables of the current run:
    %runvar contains path and file names that depend on users and input data
    
    %1) store start msg in runvar:      

    runvar.start_msg = {'Hello, to perform the Edge Response analysis you will be asked to select the following files:';
        '1) a species abundance matrix, please refer to help manual for format. abundances may be normalised';
        '2) a list of census point coordinates, either lat / lon or x/y, in a specific format, please refer to manual for format info';
        '3) a projected map of your landscape in GeoTIFF format. The map should contain continuous variations of a single variable (e.g. Tree Cover, NDVI or LAI) representing habitat suitability, scaled from 0 to 100. Categorical maps should be transformed to fit this format, please refer to manual';
        '4) optional: a binary map in GeoTIFF format where 1 = non habitat (e.g. water, road, urban, and 0 = any natural habitat (grassland, plantation, forest..)';
        'Press OK if your input files are in the correct format and you wish to proceed, otherwise Cancel'};
   
    %2) store plugin function names (if any)

    runvar.import_plugin = []; 
    %otherwise write: runvar.import_plugin = [];
    
    %also add the theshold parameter for the distance map:
    runvar.PCT = []; %computed automatically
    %otherwise write: runvar.PCT = [];
    runvar.DEI = [];
    
    %3) create temp log file path - to be used only in case of error, overwritten otherwise with userpath    
   
    %grab current time to create log file name
    tsp = datestr(now, 'yymmdd_HHMMSS');
    %and store from whereever the software was run
    runvar.log_filename = fullfile(pwd, ['Edge_response_log_' tsp '.txt']);
    
    %4) create other fields of runvar struct (empty)
    
    runvar.software_path = [];
    runvar.folder_Data_Input = [];
    runvar.folder_Data_Matlab = [];
    runvar.folder_Data_Output = [];
    
    
    runvar.PCmap_name = []; %without extension
    runvar.PCmap_matfile = []; %full path to the PCmap matlab file
    runvar.mask_name = []; %without extension
    runvar.mask_matfile = []; %full path to the habitat mask map matlab file
    
    runvar.plot_name = []; %without extension
    runvar.plot_matfile = []; %full path to the plot matlab file
    runvar.species_name = []; %without extension
    runvar.species_matfile = []; %full path to the species matlab file
    
    runvar.Binmap_matfile = [];
    runvar.Distmap_matfile = [];
    runvar.LCmap_matfile = [];
    runvar.RespPCHD_matfile = [];
    runvar.Landscape_metrics_matfile = [];
    runvar.AbPatRef_matfile = [];
    runvar.plotrange_matfile = [];
    runvar.patternpotential_matfile = [];
    runvar.species_category_matfile = [];
    

    %5) run the sequence of operations to compute the edge response
    Edge_response_ExecSequence(runvar);
    

end
        
        
