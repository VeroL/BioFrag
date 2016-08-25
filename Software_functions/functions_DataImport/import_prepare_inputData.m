%import_prepare_inputData.m

function runvar = import_prepare_inputData(runvar, param)
    
    %1) --------------------   grab names of user selected files:
    
    %select species matrix file
    filefilter = {['*' param.csv],'Comma Separated Values files (*.csv)'}; 
    DialogTitle = [param.msgbox_Title ' : Select species abundance matrix (CSV file)'];
    openingFolder = runvar.folder_Data_Input; %data input folder. relative path ok.
    [speciesFileName,speciesFilePath,fileselected] = uigetfile(filefilter,DialogTitle, openingFolder);
    %prepare to throw error if needed:
    var_to_assert = fileselected;
    erreason = 'NO SPECIES MATRIX FILE selected';
    errID = 'inputs';
    check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if no file selected
    
    %select plots file
    filefilter = {['*' param.csv],'Comma Separated Values files (*.csv)'}; 
    DialogTitle = [param.msgbox_Title ' : Select list of census point coordinates (CSV file)'];
    openingFolder = runvar.folder_Data_Input; %data input folder. relative path ok.
    [plotFileName,plotFilePath,fileselected] = uigetfile(filefilter,DialogTitle, openingFolder);
    %prepare to throw error if needed:
    var_to_assert = fileselected;
    erreason = 'NO CENSUS POINTS FILE selected';
    errID = 'inputs';
    check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if no file selected
        
    %select map tif file
    filefilter = {['*' param.tif],'GeoTIFF files (*.tif)'}; 
    DialogTitle = [param.msgbox_Title ' : Select projected land cover map (GeoTIFF file)'];
    openingFolder = runvar.folder_Data_Input; %data input folder. relative path ok.
    [tifFileName,tifFilePath,fileselected] = uigetfile(filefilter,DialogTitle, openingFolder);
    %prepare to throw error if needed:
    var_to_assert = fileselected;
    erreason = 'NO MAP GEOTIFF FILE selected';
    errID = 'inputs';
    check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if no file selected
    
    %select mask tif file if user has one
    qstring = 'Do you have a non-habitat mask map of your landscape? (habitat = 0 | non-habitat = 1, same projection and extent as landcover map)';
    yesbutton = 'Yes, select';
    nobutton = 'No, proceed without';
    button = questdlg(qstring,param.msgbox_Title,yesbutton,nobutton,yesbutton);
    switch button
        case yesbutton %user has a mask
            filefilter = {['*' param.tif],'GeoTIFF files (*.tif)'}; 
            DialogTitle = [param.msgbox_Title ' : Select projected binary mask (GeoTIFF file)'];
            openingFolder = runvar.folder_Data_Input; %data input folder. 
            [masktifFileName,masktifFilePath,fileselected] = uigetfile(filefilter,DialogTitle, openingFolder);
            %prepare to throw error if needed:
            var_to_assert = fileselected;
            erreason = 'NO MASK GEOTIFF FILE selected';
            errID = 'inputs';
            check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if no file selected
        case nobutton
            masktifFileName = [];
        otherwise 
            masktifFileName = [];
    end
    
  
    
    %2) --------------------   import geotiff map and save
    tifFileFullPath = fullfile(tifFilePath,tifFileName);
    [~,runvar.PCmap_name] = fileparts(tifFileName); %remove tif extension
    PCmap_matname =  [runvar.PCmap_name param.PCmat]; %create PCmap matlab file name
    runvar.PCmap_matfile = fullfile(runvar.folder_Data_Matlab, PCmap_matname); %create path to matlab PC map
        
    if ~exist(runvar.PCmap_matfile,'file') %if the tif file has not been processed before:
        %import the projected tif:
        [PC.map, PC.R, PC.ginfo, PC.res, ismap_geotif] = importRaster(tifFileFullPath);
        
        %check if file contains geographic info
        %prepare to throw error if needed:
        var_to_assert = ismap_geotif;
        erreason = ['The GTModelTypeGeoKey in the GeoKeyDirectoryTag is not defined in the GeoTIFF file selected: ' strrep(tifFileFullPath, '\', '/')];
        errID = 'inputs';
        check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if map is not a correct geotiff file
        
        
        %save PC structure:
        save(runvar.PCmap_matfile,'PC');
        msg = ['Imported and stored Point Cover map from file: ' tifFileFullPath]; dispwrite_log(runvar, param, msg)
    else
        msg = ['previously imported: Point Cover map tif file ' tifFileFullPath]; dispwrite_log(runvar, param, msg)
    end
    
    
    %if user has selected a mask, import it and check that the extent corresponds to the PCmap
    if ~isempty(masktifFileName)
        
        masktifFileFullPath = fullfile(masktifFilePath,masktifFileName);
        [~,runvar.mask_name] = fileparts(masktifFileName); %remove tif extension
        mask_matname =  [runvar.mask_name param.mask]; %create mask matlab file name
        runvar.mask_matfile = fullfile(runvar.folder_Data_Matlab, mask_matname); %create path to matlab mask map

        if ~exist(runvar.mask_matfile,'file') %if the tif file has not been processed before:
            %import the projected tif:
            [mask.map, mask.R, mask.ginfo, mask.res, ismap_geotif] = importRaster(masktifFileFullPath);

            %check if file contains geographic info
            %prepare to throw error if needed:
            var_to_assert = ismap_geotif;
            erreason = ['The GTModelTypeGeoKey in the GeoKeyDirectoryTag is not defined in the GeoTIFF file selected: ' strrep(masktifFileFullPath, '\', '/')];
            errID = 'inputs';
            check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if map is not a correct geotiff file
            
            %then check that the extent of the mask corresponds to the extent of the PCmap
            % load PC map (needed if was previsouly imported):
            load(runvar.PCmap_matfile);
            var_to_assert = isequal(PC.R,mask.R);
            %prepare to throw error if needed:
            erreason = 'The extent of the mask provided does not correspond to the extent of the landcover map';
            errID = 'inputs';
            check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if mask extent does not correspond to PCmap extent
            
            %then check that the mask is binary
            %note, for the biofrag DB the files provided by Marion have the following format : 0 = no data, 1 = no water, 2 = water.
            %             watermask(watermask == 0) = 2; %make the no data as water, so there are not counted as land. 
            %             watermask = watermask - 1; % so 1 means water, 0 means not water
            nbpix_inmask = numel(mask.map);
            nbpix_0n1 = sum(sum(mask.map==0)) + sum(sum(mask.map==1));
            var_to_assert = isequal(nbpix_inmask,nbpix_0n1);
            %prepare to throw error if needed:
            erreason = 'The mask provided contains other values than 0 and 1 - not binary';
            errID = 'inputs';
            check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if mask extent does not correspond to PCmap extent

            %save mask structure:
            save(runvar.mask_matfile,'mask');
            msg = ['Imported and stored mask map from file: ' masktifFileFullPath]; dispwrite_log(runvar, param, msg)
        else
            msg = ['previously imported: mask map tif file ' masktifFileFullPath]; dispwrite_log(runvar, param, msg)
        end
        
    end
    
    
    
    %3) --------------------   import plots 
    plotcsvFileFullPath = fullfile(plotFilePath,plotFileName);
    [~,runvar.plot_name] = fileparts(plotFileName); %remove csv extension
    plot_matname =  [runvar.plot_name param.plotLoc]; %create plotLoc matlab file name
    runvar.plot_matfile = fullfile(runvar.folder_Data_Matlab, plot_matname); %create path to matlab plotLoc
    
    if ~exist(runvar.plot_matfile,'file') %if the plot file has not been processed before:
        % import the CSV file, make sure plot names are string, order by natural order
        % store plots in a structure array where each array element is a plot, and fields are name, lat, lon
        [plotLoc, headerline_ok] = importPlotCoords(plotcsvFileFullPath);
        %check if import worked, only 1 possible error: can't find string names:
        %prepare to throw error if needed:
        var_to_assert = headerline_ok;
        erreason = ['cannot find strings lat and lon or x and y in headerline of file: ' strrep(plotcsvFileFullPath, '\', '/')];
        errID = 'inputs';
        check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if incorrect headerline in input plot file
   
        %if run carries on display ok
        msg = ['Imported census points from file: ' plotcsvFileFullPath]; dispwrite_log(runvar, param, msg)
    else
        load(runvar.plot_matfile);
        msg = ['previously imported: census points csv file ' plotcsvFileFullPath]; dispwrite_log(runvar, param, msg)
    end
    
    
    %4) --------------------   import species
    speciescsvFileFullPath = fullfile(speciesFilePath,speciesFileName);
    [~,runvar.species_name] = fileparts(speciesFileName); %remove csv extension
    species_matname =  [runvar.species_name param.species]; %create species matlab file name
    runvar.species_matfile = fullfile(runvar.folder_Data_Matlab, species_matname); %create path to matlab species
    
    if ~exist(runvar.species_matfile,'file') %if the species file has not been processed before:
        % import the CSV file, make sure plot names are string, order by natural order
        % store species matrix, plot names, species names, boolean of abundant species
        species = importSpeciesMatrix(speciescsvFileFullPath);
        msg = ['Imported species matrix from file: ' speciescsvFileFullPath]; dispwrite_log(runvar, param, msg)
    else
        load(runvar.species_matfile);
        msg = ['previously imported: species matrix csv file ' speciescsvFileFullPath]; dispwrite_log(runvar, param, msg)
    end
    
    
    % 5-->7) --------------------   double check species and plot files, project plots if needed, save

    if ~exist(runvar.plot_matfile,'file') || ~exist(runvar.species_matfile,'file') %if plot or species or both not stored as mat files:
        
        %5) --------------------   check plot names correspond in plot and species files
        samenb_plots = isequal(length(plotLoc),length(species.plotNames));
        %prepare to throw error if needed:
        var_to_assert = samenb_plots;
        erreason = ['mismatch in number of census points in file ' strrep(plotcsvFileFullPath, '\', '/') ' and in file ' strrep(speciescsvFileFullPath, '\', '/')];
        errID = 'inputs';
        check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if plot names do not correspond
        
        plotnamesok = all(strcmp({plotLoc.name}',species.plotNames)) == 1;
        %prepare to throw error if needed:
        var_to_assert = plotnamesok;
        erreason = ['mismatch in census point names in file ' strrep(plotcsvFileFullPath, '\', '/') ' and in file ' strrep(speciescsvFileFullPath, '\', '/')];
        errID = 'inputs';
        check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if plot names do not correspond
        
    
        %6) --------------------   project plots
        % load PC map to project plot coord on it (get row and col):
        load(runvar.PCmap_matfile);
        %project plots if needed (if lat|lon was given), then compute row and col of plots, and add coordinates into plotLoc structure
        [plotLoc, mapcoord_match] = projectPlotsLatLon(plotLoc, PC);  %#ok<ASGLU>
        %if not all plots fall on map provided send error:
        %prepare to throw error if needed:
        var_to_assert = mapcoord_match;
        erreason = ['not all census points coordinate fall within the extent of the map provided: ' strrep(tifFileFullPath, '\', '/')];
        errID = 'inputs';
        check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if plots fall outside the map
      

        %7) --------------------   save plot and species as matlab files
        %save plotLoc
        save(runvar.plot_matfile,'plotLoc'); 
        %save species matrix
        save(runvar.species_matfile,'species'); 
        msg = ['Files ' plotcsvFileFullPath ' and ' speciescsvFileFullPath ' have been successfully checked and stored']; dispwrite_log(runvar, param, msg)
    end
    
        
    msg='* Edge response is computed using the following files:'; dispwrite_log(runvar, param, msg)
    msg=['-> Point Cover map: ' tifFileFullPath]; dispwrite_log(runvar, param, msg)
    if ~isempty(masktifFileName)
        msg=['-> mask map: ' masktifFileName]; dispwrite_log(runvar, param, msg)
    end
    msg=['-> Census point coordinates: ' plotcsvFileFullPath]; dispwrite_log(runvar, param, msg)
    msg=['-> Species abudance matrix: ' speciescsvFileFullPath]; dispwrite_log(runvar, param, msg)
 

end