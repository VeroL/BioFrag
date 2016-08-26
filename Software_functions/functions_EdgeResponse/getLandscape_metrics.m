% getLandscape_metrics

%compute the total habitat amout in landscape (habitat defined by PC >60, soft threshold), in hectares
%habitat configuration index (total edge effects but rescaled)
%total edge effects on habitat

%stores in runvar.Landscape_metrics_matfile as a structure with 3 fields : PC_amount, PC_CI, PC_Corescore

function runvar = getLandscape_metrics(runvar, param)

    %create matlab filename for landscape metrics
    runvar.Landscape_metrics_matfile = fullfile(runvar.folder_Data_Matlab, [runvar.PCmap_name '_DEI_' num2str(runvar.DEI) param.Landscape_metrics]);
    % to store Landscape_metrics, structure (1) with fields PC_amount, PC_CI, PC_Corescore
    
    if ~exist(runvar.Landscape_metrics_matfile,'file')
        
        
        %load PC map, EI map
        load(runvar.LCmap_matfile);%for HD map
        load(runvar.PCmap_matfile); % for FI
        
        %assume no water / no non-habitat
        mask.map = zeros(size(PC.map)); 
        %but overwrite / load the mask if it exists 
        if ~isempty(runvar.mask_matfile)
            if exist(runvar.mask_matfile,'file')
                clear mask
                load(runvar.mask_matfile);
            end
        end
        
        %variable extraction:
        PCmap = PC.map;
        HDmap = HD(1).map;        
        res = PC.res; %resolution
        
        %1) ------------------- compute Habitat Amount in ha
        
        %make sigmoid function to account for denser pixels more, saturates arount 80, (so pixels with low TC count for less area)
        infpointPCf = 60;  %soft threshold parameters for TC
        ratePC = 0.25;
        
        PCmap_thresh = 1./(1+exp(-ratePC*(PCmap - infpointPCf))); 

        PC_amount = sum(PCmap_thresh(:) .* ~mask.map(:)); %in pixels * density, i.e. nb of trees,  makes sure water = 0
        Landscape_metrics.PC_amount = PC_amount .* res^2 / 10000; %in hectares 
        
        
        %2) ------------------- compute Configuration Index
        
        PCmap = PCmap .* ~mask.map; % make sure water pixels have a score of 0
        
        %---Clip
        binmap = double(PCmap > 30); %low thresold, just make sure theres not lots of matrix around
        %find habitat extreme positions to clip the map around:
        [rowhab, colhab] = find(binmap);        
        %clip map using min and max positions of habitat, no padding :        
        PCmapclip = PCmap(min(rowhab):max(rowhab),min(colhab):max(colhab));
        size_PCmapclip = numel(PCmapclip); %number of pixels
        %---Compute EImap with DEI that depends on image size:
        DEI_forCI = sqrt(size_PCmapclip)/10 * res; % the / 10 corresponds to a reference image of 10km side
        EImap_forCI = MakeEIMap(PCmapclip, res, DEI_forCI);        
        %---Compute the score with clip map and variable EImap_forCI:
        Landscape_metrics.PC_CI = get_habitat_score(PCmapclip, EImap_forCI, infpointPCf, ratePC);
        
        
        %3 ------------------- compute Core Score
        
        Landscape_metrics.PC_Corescore = get_habitat_score(PCmap, HDmap, infpointPCf, ratePC);
        
        
        %-------------------store .mats
        save(runvar.Landscape_metrics_matfile ,'Landscape_metrics');

        msg=['computed Landscape metrics for ' runvar.PCmap_name]; dispwrite_log(runvar, param, msg)

    else
        msg=['previously computed: Landscape metrics ' runvar.Landscape_metrics_matfile]; dispwrite_log(runvar, param, msg)

    end %end of condition on output existence

end %end function



% function to compute the core to habitat ratio (core score)
function  PC_Corescore = get_habitat_score(PCmap, EImap, infpointPCf, ratePC)
 
    %soft threshold parameters for PC given in input
    %soft threshold parameters for EI:
    infpointEIc = 20;
    rateEI = 0.35;
                
    %forest score map (soft threshold forest)
    HabitatScoremap = 1./(1+exp(-ratePC*(PCmap - infpointPCf))); %between 0 and 1
    
    %core score map  (soft threshold EI)
    Corescoremap = 1./(1+exp(rateEI*(abs(EImap) - infpointEIc ))); %between 0 and 1
    
    %forest core score map
    effectiveHabitatScoremap = Corescoremap .* HabitatScoremap;  %between 0 and 1 

    PC_Corescore = sum(effectiveHabitatScoremap(:))/sum(HabitatScoremap(:));   
    
end



    
    