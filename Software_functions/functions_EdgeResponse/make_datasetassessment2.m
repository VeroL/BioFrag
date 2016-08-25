% dataset assessment :

%compute plot range

%read in plot range

%read in data summary like in report csv

%read in watermask
%do area calculations  (now in compute areas)

% P(species category  is correct | dataset suitable is suitbale)

% P(Sc) = P(Sc AND Ds) + P(Sc AND ¬Ds) 
% P(Sc) = P(Sc | Ds)*P(Ds) + 0;

function runvar = make_datasetassessment2(runvar, param)

    %create filenames for plotrange and patternpotential
    runvar.plotrange_matfile = fullfile(runvar.folder_Data_Matlab, [runvar.plot_name '_DEI_' num2str(runvar.DEI) param.plotrange]);
    runvar.patternpotential_matfile = fullfile(runvar.folder_Data_Matlab, [runvar.plot_name '_DEI_' num2str(runvar.DEI) param.patternpotential]);
    
    
    if ~exist(runvar.plotrange_matfile,'file') || ~exist(runvar.patternpotential_matfile,'file') 
        
        
        %########################## PLOT RANGE ######################################

        %INPUT VALUES--------------------------------------------------------------

        %min valid range of plots on PCHD graph to assess which patterns can be evaluated:
        minplotrange_HDlarge = 40;
        minplotrange_HD = 3;
        minplotrange_PC = 60;

        %penalty for dataset rating:
        penalty = 9/10;
        penaltyDEI = 7/10;
        
        %list of categories we deem essential. Plots should be distributed so that species in these categories can be assessed.
        essential_subcats = {   %modified version of subcat, without habitat tag. when to comparing to pattern's subcat the tag is stripped (in pothabs_filt)                                   
        'Forest Core 0'
        'Matrix Core 0'
        'Generalist Core 0'

        'Forest Edge T3 highside'
        'Matrix Edge T3 lowside'
        'Generalist Edge T3'

        'Forest Edge T2 highside'
        'Matrix Edge T2 lowside'
        'Generalist Edge T2'

        'Forest noPref inf'
        'Matrix noPref inf'
        'Generalist noPref'
        };

        %------------------------load species smooth abundance PCHD and local cover, and pattern training set
        load(runvar.RespPCHD_matfile);
        load(runvar.LCmap_matfile);
        load(runvar.AbPatRef_matfile);

        %get index of optimal DEI:
        DEIopt = RespPCHD(1).DEIopt;


        % ---------------- find set of reachable patterns for dataset (depends on plot distribution and landscape (for non opt DEI)) -----------
        %---------------------------------------------------------------------------------------------------------------------------------------

        % do this for all DEIs
        % stats on reachable patterns in dataset

        %variable memory allocation
        nbplots = length(RespPCHD(1).plot_PointCover);
        nbpatterns = length(AbPatRef);
        nbDEIs = length(RespPCHD);

        pattern_ispotential = false(nbpatterns,nbDEIs);
        plot_patternAb = zeros(nbplots,nbpatterns,nbDEIs); %containg abundance values for all plots, for all patterns and DEI

        plot_range(1).propShapes = 0;
        plot_range(1).hasForestCore = 0;
        plot_range(1).hasForestEdge = 0; 
        plot_range(1).hasMatrixCore = 0;
        plot_range(1).hasMatrixEdge = 0; 
        plot_range(1).DEIopt = 0;
        plot_range(1).rating = 0;
        plot_range(1).DEI = 0;
        plot_range(2:nbDEIs) = plot_range(1);


        for indDEI = 1 : nbDEIs %loop on DEI values

            %get plot ranges in PC and HD
            PCrange = max(RespPCHD(indDEI).plot_PointCover) - min(RespPCHD(indDEI).plot_PointCover);
            HDrange = max(RespPCHD(indDEI).plot_HabitatDisturbance) - min(RespPCHD(indDEI).plot_HabitatDisturbance);

            %get plot indices in PCHD graph , 201 rows, 101 cols
            %requires round values between 1 and 201 for HD, 1 and 101 for PC
            %HD is real between -100 and 100, PC is real between 0 and 100
            [graph_nrow, graph_ncol] = size(AbPatRef(1).matrix);
            indPlotHD = round(RespPCHD(indDEI).plot_HabitatDisturbance)+101; %in the graph the top rows are the lowest HD values
            indPlotPC = round(RespPCHD(indDEI).plot_PointCover)+1;
            indPlot_ongraph = sub2ind([graph_nrow, graph_ncol], indPlotHD , indPlotPC); %line indices


            for indPat = 1 : nbpatterns %loop on abundance patterns

                %get zone (max ab or not) of each plot
                plot_zone = AbPatRef(indPat).matrix_maxbin(indPlot_ongraph);

                %get plot HabitatDisturbance values and HDrange in each zone:
                % 0 (background), 1 (high abundance) and 2 (second high abundance zone for generalists)
                PCrange_inzone = zeros(1,3);
                HDrange_inzone = zeros(1,3);
                for z = 1:3
                    %PC range in zone:
                    PCplot_inzone = RespPCHD(indDEI).plot_PointCover(plot_zone == z-1);
                    PCrange_inzone_tmp = max(PCplot_inzone) - min(PCplot_inzone);
                    if ~isempty(PCrange_inzone_tmp)
                        PCrange_inzone(z) = PCrange_inzone_tmp;
                    end

                    %HD range in zone:
                    HDplot_inzone = RespPCHD(indDEI).plot_HabitatDisturbance(plot_zone == z-1);
                    HDrange_inzone_tmp = max(HDplot_inzone) - min(HDplot_inzone);
                    if ~isempty(HDrange_inzone_tmp)
                        HDrange_inzone(z) = HDrange_inzone_tmp;
                    end
                end

                %check potentiality of pattern depending on category
                switch AbPatRef(indPat).Category
                    case {'Forest Core','Forest Edge','Matrix Core','Matrix Edge'}
                        %must have HD range over 5 in high zone (2) and outside (1)
                        if HDrange_inzone(1) > minplotrange_HD && HDrange_inzone(2) > minplotrange_HD, pattern_ispotential(indPat,indDEI) = 1; end
                    case {'Forest noPref','Matrix noPref'}
                        %must have HD range > 40 in high zone 
                        if HDrange_inzone(2) > minplotrange_HDlarge, pattern_ispotential(indPat,indDEI) = 1; end
                    case {'Generalist Core'}
                        if max(AbPatRef(indPat).matrix_maxbin(:)) == 2  % if patterns has 2 max zones
                            %must have HD range over 5 in both high zones and outside
                            if HDrange_inzone(1) > minplotrange_HD && HDrange_inzone(2) > minplotrange_HD && HDrange_inzone(3) > minplotrange_HD, pattern_ispotential(indPat,indDEI) = 1; end
                        else                                            % if patterns has 1 max zone
                            %must have HD range over 5 in high zone and outside AND PCrange_inzone over 50 in the one high zone and outside
                            if PCrange_inzone(1) > minplotrange_PC && PCrange_inzone(2) > minplotrange_PC && HDrange_inzone(1) > minplotrange_HD && HDrange_inzone(2) > minplotrange_HD
                                pattern_ispotential(indPat,indDEI) = 1; 
                            end
                        end
                    case {'Generalist Edge'}
                        %must have HD range over 5 in both high zones and outside AND PCrange over 50
                        if PCrange > minplotrange_PC && HDrange_inzone(1) > minplotrange_HD && HDrange_inzone(2) > minplotrange_HD && HDrange_inzone(3) > minplotrange_HD, pattern_ispotential(indPat,indDEI) = 1; end
                    case 'Generalist noPref'
                        %must have global HD range over 60, and PC range over 50
                        if PCrange > minplotrange_PC && HDrange > 2*minplotrange_HDlarge, pattern_ispotential(indPat,indDEI) = 1; end
                    case {'Nonabundant','Unknown','Absent'}
                        %unknown and Nonabundant are always a potential patterns
                        pattern_ispotential(indPat,indDEI) = 1; 
                    otherwise
                        warning(['cannot find pattern category: ' AbPatRef(indPat).Category]);
                end

                % % -------------------------------------------------------image verification of code
                %         if indDEI==2  %&& (strcmp(AbPatRef(indPat).Category,'Forest Edge') || strcmp(AbPatRef(indPat).Category,'Generalist Edge'))
                %             if AbPatRef(indPat).isScalemax
                %                 %if strcmp(AbPatRef(indPat).Category,'Generalist Core')
                %                 %figure;
                %                 imagesc(AbPatRef(indPat).matrix_maxbin);
                %                 if pattern_ispotential(indPat,indDEI)
                %                     colormap([1 1 1; 0 1 0; 0 0.5 0]);
                %                 else
                %                     colormap([1 1 1; 1 0 0; 0.5 0 0]);
                %                 end
                %                 hold on
                %                 plot(indPlotPC, indPlotHD, 'bx'); 
                %                 plot([1 101],[101 101],'k'); plot([1 101],[101 51],'k--'); plot([1 101],[151 101],'k--');
                %                  plot([1 101],[101 1],'k-'); plot([1 101],[201 101],'k-');
                %                 hold off
                %                 set(gca,'YDir','normal'); axis([1 101 1 201]); axis image
                %                 title([AbPatRef(indPat).ScaleCategory ' - is valid ? ' num2str(pattern_ispotential(indPat,indDEI)) ', DEI =' num2str(RespPCHD(indDEI).DEI)]);
                %                 %LCrange_inzone
                %                 waitforbuttonpress;
                %                 %end
                %             end
                %         end
                % % ---------------------------------------------------------------------------------------------------

                
            end %end of first loop on patterns --------------------------------------------------------------

            %---------------- second assessment of potential patterns
        
            % Forest noPref patterns -> need forest core + forest edge 
            % Matrix noPref patterns -> need Matrix core + Matrix edge 
            % Generalist core -> need matrix core + forest core
            % generalist edge -> need matrix edge + forest edge
            % generalsit nopref -> need matrix core matrix edge forest core forest edge (i.e. need gen core and gen edge)

            %after first selection compute the list of reachable main categories:
            reachable_main_categories = unique({AbPatRef(pattern_ispotential(:,indDEI)).Category}');

            % Forest noPref patterns -> need forest core + forest edge  | if not both reachbale set forest noPref to not potential
            if ~ ( ismember({'Forest Core'},reachable_main_categories) && ismember({'Forest Edge'},reachable_main_categories) )
                pattern_ispotential(strcmp({AbPatRef.Category}','Forest noPref'),indDEI) = 0;
            end
            % Matrix noPref patterns -> need Matrix core + Matrix edge  | if not both reachbale set Matrix noPref to not potential
            if ~ ( ismember({'Matrix Core'},reachable_main_categories) && ismember({'Matrix Edge'},reachable_main_categories) )
                pattern_ispotential(strcmp({AbPatRef.Category}','Matrix noPref'),indDEI) = 0;
            end
            % Generalist core -> need matrix core + forest core  | if not both reachbale set Generalist Core to not potential
            if ~ ( ismember({'Forest Core'},reachable_main_categories) && ismember({'Matrix Core'},reachable_main_categories) )
                pattern_ispotential(strcmp({AbPatRef.Category}','Generalist Core'),indDEI) = 0;
            end
            % generalist edge -> need matrix edge + forest edge  | if not both reachbale set Generalist Edge to not potential
            if ~ ( ismember({'Forest Edge'},reachable_main_categories) && ismember({'Matrix Edge'},reachable_main_categories) )
                pattern_ispotential(strcmp({AbPatRef.Category}','Generalist Edge'),indDEI) = 0;
            end

            % generalsit nopref -> need Generalist core + generalist edge  | if not both reachbale set generalsit nopref to not potential
            if ~ ( ismember({'Forest Core'},reachable_main_categories) && ismember({'Matrix Core'},reachable_main_categories) && ...
                    ismember({'Forest Edge'},reachable_main_categories) && ismember({'Matrix Edge'},reachable_main_categories) )
                pattern_ispotential(strcmp({AbPatRef.Category}','Generalist noPref'),indDEI) = 0;
            end

            %recompute main categories
            reachable_main_categories = unique({AbPatRef(pattern_ispotential(:,indDEI)).Category}');
            
            %---------------------------------------new loop on pattern after pattern_ispotential has benn corrected: 
            for indPat = 1 : nbpatterns %loop on abundance patterns
                %if pattern is potential, save its abundance values for each plot:
                if pattern_ispotential(indPat,indDEI)
                    plot_patternAb(:,indPat,indDEI) = AbPatRef(indPat).matrix(indPlot_ongraph);
                end
            end
            %---------------------------------------           
            
            %--------------------------------dataset stats on plot distribution:

            %get list of potential pattern habitat: forest matrix or generalist (remove dense/ sparse tag)
            pothabs = {AbPatRef(pattern_ispotential(:,indDEI)).Habitat}';
            pothabs_filt = cellfun(@(x) strtok(x, ' '),pothabs,'UniformOutput',false); %take first token, i.e. remove the habitat tag (e.g. 'Sp')

            %get list of potential pattern EI preference (with T" highside tags, etc)
            poteips = {AbPatRef(pattern_ispotential(:,indDEI)).EIPref}';

            %get list of reachable pattern subcategories by associating habitat and EIpref pattern lists
            reachable_subcats = cellfun(@(x,y) [x ' ' y],pothabs_filt,poteips,'UniformOutput',false);
            %list of unique reachable subcategories
            reachable_subcats = unique(reachable_subcats);
            %compare with list of essential categories (intersection)
            reachable_essential_subcats = intersect(essential_subcats,reachable_subcats);
            %compute the ratio of reachable essential cateogirs to the number of essential categories:
            propshapes = length(reachable_essential_subcats) / length(essential_subcats);

            %store:
            plot_range(indDEI).propShapes = propshapes;            
            
            
            

            %booleans, whether plots have been taken in forest core and edge, and matrix core and edge       
            plot_range(indDEI).hasForestCore = any(strcmp('Forest Core',reachable_main_categories));
            plot_range(indDEI).hasForestEdge = any(strcmp('Forest Edge',reachable_main_categories)); 
            plot_range(indDEI).hasMatrixCore = any(strcmp('Matrix Core',reachable_main_categories));
            plot_range(indDEI).hasMatrixEdge = any(strcmp('Matrix Edge',reachable_main_categories)); 

            plot_range(indDEI).DEIopt = DEIopt;

            plot_range(indDEI).rating = plot_range(indDEI).propShapes ...
                * (penalty*(~plot_range(indDEI).hasForestCore) + plot_range(indDEI).hasForestCore) ...
                * (penalty*(~plot_range(indDEI).hasForestEdge) + plot_range(indDEI).hasForestEdge) ...
                * (penalty*(~plot_range(indDEI).hasMatrixCore) + plot_range(indDEI).hasMatrixCore) ...
                * (penalty*(~plot_range(indDEI).hasMatrixEdge) + plot_range(indDEI).hasMatrixEdge) ...
                * (penaltyDEI*(DEIopt < 500) + (DEIopt >= 500));



           %store DEI for which assessment has been computed
           plot_range(indDEI).DEI = RespPCHD(indDEI).DEI;
           
           
                        

        end %end of DEI loop ------------------------------------



        %store .mats
        save( runvar.plotrange_matfile,'plot_range');
        save( runvar.patternpotential_matfile,'pattern_ispotential','plot_patternAb');


        msg=['computed plot ranges, pattern potential and abvalues saved as mat for ' runvar.plot_name]; dispwrite_log(runvar, param, msg)

        %     pattern_ispotential = false(nbpatterns,nbDEIs); %boolean, wether a pattern (row) is potential for this dataset (file) and this dei (col)
        %     plot_patternAb = zeros(nbplots,nbpatterns,nbDEIs); %containg abundance values for all plots, for all patterns and DEI  
        
    else
        msg=['previously computed: plot ranges, pattern potential ' runvar.plotrange_matfile]; dispwrite_log(runvar, param, msg)

    end %end of condition on output existence
    
end


   
    
    
    
    
    
    
    
    
    
    