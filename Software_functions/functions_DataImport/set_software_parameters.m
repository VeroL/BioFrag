%set_software_parameters

%param structure contains naming conventions of files and folder and messages
%these do not vary

function param = set_software_parameters

    %software name for title of message windows:
    param.msgbox_Title = 'Edge response';
    
    %error msgs:
    param.runaborted = ['!!! Error: ' param.msgbox_Title ' software cannot complete because: '];
    param.errorID = 'EdgeResponseSoft:';
    
    
    %folder names:
    param.foldername_Data_Input = 'Edge_response_Data_Input';
    param.foldername_Data_Matlab = 'Edge_response_Data_Matlab';
    param.foldername_Data_Output = 'Edge_response_Data_Output';
   
    %list filenames conventions:
    
    %input file names:
    param.tif = '.tif';
    param.masktif = '_mask.tif';
    param.csv = '.csv';    
        
    %matlab working file names:
    param.PCmat = '_PointCover.mat';
    param.mask = '_mask.mat';
    param.plotLoc = '_plots.mat';
    param.species = '_species.mat';
    param.PCBin = '_PointCoverBin.mat';
    param.Distance2nE = '_Distance2nE.mat';
    param.LocalCover = '_LocalCover.mat';
    param.RespPCHD = '_RespPCHD.mat';
    param.AbPatRef = 'AbPatRef.mat';
    param.plotrange='_plotrange.mat';    
    param.patternpotential='_patternpotential.mat';
    param.species_category='_species_category.mat';
        
    %output file names:
    param.PCBintif = '_binary.tif';
    param.Distance2nEtif = '_distance2nearestEdge.tif';
    param.EImaptif = '_EImap.tif';
    param.censuspoints_propertiesCSV = '_censuspoints_properties.csv';
    param.smoothAbundancePCEICSV = '_smoothAbundancePCEI.csv';
    param.species_categoryCSV = '_species_category.csv';
    param.species_category_PP_CSV = '_species_category_posteriorprob.csv';  
    
   
end