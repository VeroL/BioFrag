%generate software file outputs: 

%geotiff of Binary map
%geotiff of distance map
%geotiff of EI map

%csv containing census plots properties: plotName x y row col distance2nearestEdge PC EI
%csv containing species matrix smoothed with respect to PC and EI at census locations
%csv containing list of species with the following attributes: 'Species name','Habitat','EIPref','Category','Posterior prob', 'Pattern Label', 'DEI opt','mean abundance', ...
%             'median abundance','fragmentation impact','EI sensitivity','dataset rating'
%csv contanining list of species and their probablity to belong to each of the main categories

%tif file of EI, csv with categories, FI and dataset rating, csv with distance

function generate_outputs(runvar, param)

    %input values
    overwriteL = 0; %1 to overwrite tif maps
    overwriteD = 0; %1 to overwrite csv records
    
    %count number of files generated:
    nbgenfiles = 0;

    %load data
    load(runvar.species_category_matfile); %for species response
    load(runvar.species_matfile); %for species names
    load(runvar.LCmap_matfile); %for EI map
    load(runvar.plotrange_matfile); %for dataset rating
    load(runvar.Binmap_matfile); %for binary map
    load(runvar.Distmap_matfile); %distance map
    load(runvar.PCmap_matfile); %to get tif info
    load(runvar.AbPatRef_matfile); %to get species category list
    main_category_list = AbPatRef_properties.main_category_list;
    load(runvar.plot_matfile); %for plot coordinates and name + store distances + plot row and cols
    load(runvar.RespPCHD_matfile); %for plot PC and HD and smooth abs 
    load(runvar.Landscape_metrics_matfile); %for landscape metrics
    
    nbspecs = length(species.names);
    
    %grab index of opt DEI
    indDEI = SpeciesEEresp(1).DEIoptind; %i.e. only using optimal DEI
    DEI = HD(indDEI).DEI;
    
    
    %------------------------------------Binary map tif :
    
    Binarytif = fullfile(runvar.folder_Data_Output, [runvar.PCmap_name param.PCBintif]);
    
    %if the file does not exist - or if it does but overwrite is set to 1
    if ~exist(Binarytif,'file') || overwriteL
    
        R = PC.R;
        ginfo = PC.ginfo;

        writeRastertif(Binarytif,PCBin.map,R,ginfo);

        msg=['Binary map exported to tif for ' runvar.PCmap_name]; dispwrite_log(runvar, param, msg)
        nbgenfiles = nbgenfiles + 1;
    end
    
    
    %------------------------------------Distance map tif :
    
    Distancetif = fullfile(runvar.folder_Data_Output, [runvar.PCmap_name param.Distance2nEtif]);
    
    %if the file does not exist or if it does but overwrite is set to 1
    if ~exist(Distancetif,'file') || overwriteL
    
        R = PC.R;
        ginfo = PC.ginfo;

        writeRastertif(Distancetif,D2nE.map,R,ginfo);

        msg=['Distance map exported to tif for ' runvar.PCmap_name]; dispwrite_log(runvar, param, msg)
        nbgenfiles = nbgenfiles + 1;
    end
           
    
    %------------------------------------EI map tif for opt DEI 
    
    EImaptif = fullfile(runvar.folder_Data_Output, [runvar.PCmap_name '_DEI' num2str(DEI) param.EImaptif]);
    
    %if the file does not exist or if it does but overwrite is set to 1
    if ~exist(EImaptif,'file') || overwriteL
    
        R = PC.R;
        ginfo = PC.ginfo;
        
        EImap = HD(indDEI).map;

        writeRastertif(EImaptif,EImap,R,ginfo);

        msg=['EI map exported to tif for ' runvar.PCmap_name]; dispwrite_log(runvar, param, msg)
        nbgenfiles = nbgenfiles + 1;
    end
    
    
    
    
    %------------------------------------ plot properties CSV (+ plot PC, EI, dist, row and col, X and Y):
      
    % CSV plotfilename_DEI300_censuspoints_properties.csv ------- plot names, x, y, row, col, dist, pc, ei
    
    %filename:
    csvfile = fullfile(runvar.folder_Data_Output,[runvar.plot_name '_DEI' num2str(DEI) param.censuspoints_propertiesCSV]);
    
    %if the file does not exist - or if it does but overwrite is set to 1
    if ~exist(csvfile,'file') || overwriteD
    
        
        %------------save csv file - plot properties

        headerline = {'Plot name','X','Y','row','column','Distance','Point Cover','Edge Influence'};

        %raw matrix
        cellSummary = num2cell([ [plotLoc.X]' [plotLoc.Y]' ...
            [plotLoc.pR]' [plotLoc.pC]' ...
            [plotLoc.dist]' [RespPCHD(indDEI).plot_PointCover] [RespPCHD(indDEI).plot_HabitatDisturbance] ]);
 
        %add plot name column
        cellSummary = [{plotLoc.name}' cellSummary]; %#ok<*AGROW>
        %add header lines:
        cellSummary = [headerline; cellSummary];
        
        %export as CSV file
        exportas_CSV(cellSummary,csvfile);
        clear cellSummary;

        msg=['census points properties saved as CSV for ' runvar.plot_name]; dispwrite_log(runvar, param, msg)
        nbgenfiles = nbgenfiles + 1;
    end
    
    
    
    
    %------------------------------------ smooth species response wrt PCHD CSV  (+ plot names):
    
    % CSV speciesfilename_DEI300_smoothAbundancePCEI.csv ------- smoothed species matrix and plot names
    
    %filename:
    csvfile = fullfile(runvar.folder_Data_Output,[runvar.species_name '_DEI' num2str(DEI) param.smoothAbundancePCEICSV]);
    
    %if the file does not exist - or if it does but overwrite is set to 1
    if ~exist(csvfile,'file') || overwriteD    
        
        %------------save csv file - Resp to PCHD 

        headerline = ['Plot name' species.names'];

        %raw matrix
        cellSummary = num2cell([RespPCHD(indDEI).spemat_smoothPCHD]);

        %add plot name column
        cellSummary = [{plotLoc.name}' cellSummary]; %#ok<*AGROW>
        %add header lines:
        cellSummary = [headerline; cellSummary];
        
        %export as CSV file
        exportas_CSV(cellSummary,csvfile);
        clear cellSummary;

        msg=['species abundance smoothed wrt PC and EI saved as CSV for ' runvar.species_name]; dispwrite_log(runvar, param, msg)
        nbgenfiles = nbgenfiles + 1;
    end
    
    
    
    
    %------------------------------------species category CSV
   
    % CSV speciesfilename_DEI300_species_category.csv ------- category and FI for DEI optimum +dataset rating
    
    %filename:
    csvfile = fullfile(runvar.folder_Data_Output,[runvar.species_name '_DEI' num2str(DEI) param.species_categoryCSV]);
    
    %if the file does not exist or if it does but overwrite is set to 1
    if ~exist(csvfile,'file') || overwriteD  
        
        %------------save csv file - species category 

        headerline = {'Species name','Habitat','EIPref','Category','Posterior prob', 'Pattern Label', 'DEI opt','mean abundance', ...
             'median abundance','fragmentation impact','EI sensitivity','dataset rating'} ;

        %raw matrix
        %store categories and their properties for optimum DEI

        cellSummary = num2cell([ SpeciesEEresp(indDEI).Category_prob ...
            SpeciesEEresp(indDEI).ind_AbPatRed_closestpattern ...
            repmat(DEI,nbspecs,1) ...
            SpeciesEEresp(indDEI).mean_abundance ...
            SpeciesEEresp(indDEI).median_abundance ...
            SpeciesEEresp(indDEI).Fragmentation_Impact_nowater ...
            SpeciesEEresp(indDEI).EIsensitivity ...
            repmat(plot_range.rating,nbspecs,1)]);

        %add species name column
        cellSummary = [species.names SpeciesEEresp(indDEI).Habitat SpeciesEEresp(indDEI).EIPref ...
            SpeciesEEresp(indDEI).Category cellSummary]; %#ok<*AGROW>
        %add header lines:
        cellSummary = [headerline; cellSummary];

        %export as CSV file
        exportas_CSV(cellSummary,csvfile);
        clear cellSummary;

        msg=['species category DEIopt saved as CSV for ' runvar.species_name]; dispwrite_log(runvar, param, msg)
        nbgenfiles = nbgenfiles + 1;
    end
       
    
    
    
    %------------------------------------species posterior probabilities CSV
    
    % CSV speciesfilename_DEI300_species_category_posteriorprob.csv ------- posterior prob for DEI optimum
    
    %filename:
    csvfile = fullfile(runvar.folder_Data_Output, [runvar.species_name '_DEI' num2str(DEI) param.species_category_PP_CSV]);
       
    
    %if the file does not exist or if it does but overwrite is set to 1
    if ~exist(csvfile,'file') || overwriteD               
                
        %------------save csv file - posterior probs 
        
        headerline = ['Species name' main_category_list'] ;

        %raw matrix
        cellSummary = num2cell(SpeciesEEresp(indDEI).Category_prob_matrix);

        %add species name column
        cellSummary = [species.names cellSummary]; %#ok<*AGROW>
        %add header lines:
        cellSummary = [headerline; cellSummary];

        %export as CSV file
        exportas_CSV(cellSummary,csvfile);
        clear cellSummary;

        msg=['species POSTERIOR Matrix DEIopt saved as CSV for ' runvar.species_name]; dispwrite_log(runvar, param, msg)
        nbgenfiles = nbgenfiles + 1;
    end
    
    
    
    %------------------------------------ landscape metrics and landscape input filenames :
      
    % CSV PCmapfilename_Landscape_metrics.csv ------- PCmap file name, watermask file name, DEI, habitat amount (in ha), configuration index, core score
    
    %filename:
    csvfile = fullfile(runvar.folder_Data_Output,[runvar.PCmap_name '_DEI' num2str(DEI) param.Landscape_metricsCSV]);
    
    %if the file does not exist - or if it does but overwrite is set to 1
    if ~exist(csvfile,'file') || overwriteD    
        
        %------------save csv file - landscape metrics

        headerline = {'Land Cover map','Non-habitat mask','DEI (in m)','habitat amount (in ha)','configuration index','core score'};

        %raw matrix
        cellSummary = {runvar.PCmap_name, runvar.mask_name, DEI, Landscape_metrics.PC_amount, Landscape_metrics.PC_CI, Landscape_metrics.PC_Corescore};
        
        %add header lines:
        cellSummary = [headerline; cellSummary];
        
        %export as CSV file
        exportas_CSV(cellSummary,csvfile);
        clear cellSummary;

        msg=['Landscape metrics saved as CSV for ' runvar.PCmap_name]; dispwrite_log(runvar, param, msg)
        nbgenfiles = nbgenfiles + 1;
    end
    
    
    
    
    %----------------------- final msg
    
     msg=[num2str(nbgenfiles) ' output files were generated in total, overwriteL = ' num2str(overwriteL) ' and overwriteD = ' num2str(overwriteD) '.']; dispwrite_log(runvar, param, msg)
    
    
end %end of function

                