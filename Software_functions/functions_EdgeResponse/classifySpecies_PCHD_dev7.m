% classifySpecies_PCHD_dev7

%also makes figures. we have to make figures here because don't want to store a map of predicted abundance for each species.

function runvar = classifySpecies_PCHD_dev7(runvar, param)

    %create filenames for species Edge response category (and FI)
    runvar.species_category_matfile = fullfile(runvar.folder_Data_Matlab, [runvar.species_name '_DEI_' num2str(runvar.DEI) param.species_category]);
    
    if ~exist(runvar.species_category_matfile,'file')

        %load species, smooth abundance PCHD and local cover, and pattern training set, PC map, and plot range and patternpotential
        load(runvar.species_matfile); %for species names (and raw ab to plot distance)
        load(runvar.RespPCHD_matfile); %for smooth abundance and DEI optimum
        load(runvar.LCmap_matfile);%for LC map + HD map
        load(runvar.AbPatRef_matfile); %for training set
        load(runvar.PCmap_matfile); % for FI
        load(runvar.patternpotential_matfile); %for subset of valid patterns and their abundances %loads pattern_ispotential and plot_patternAb
        
        %assume no water / no non-habitat
        mask.map = zeros(size(PC.map)); 
        %but overwrite / load the mask if it exists 
        if ~isempty(runvar.mask_matfile)
            if exist(runvar.mask_matfile,'file')
                clear mask
                load(runvar.mask_matfile);
            end
        end
        
        %to make the figure also load:        
        load(runvar.plot_matfile); %for plot coordinates
        load(runvar.plotrange_matfile); %for dataset rating

        
        % ------------------------------------------------------------------------ find best category for each species (classify)---------------
        %---------------------------------------------------------------------------------------------------------------------------------------

        %variable extraction
        %categories
        main_category_list = AbPatRef_properties.main_category_list;

        %get index of optimal DEI:
        indDEIopt = RespPCHD(1).indDEIopt;
        %DEIopt = RespPCHD(1).DEIopt;

        %PC map
        PCmap = PC.map;    

        %lengths:
        nbmaincats = length(main_category_list);
        nbspecs = length(species.names);
        nbpatterns = length(AbPatRef);
        nbDEIs = length(LC);


        %memory allocation
        SpeciesEEresp(1).Category = cell(nbspecs,1); %forest core, edge t2, t2 ....
        SpeciesEEresp(1).Habitat = cell(nbspecs,1);%forest matrix gen unknown nonabundant
        SpeciesEEresp(1).EIPref = cell(nbspecs,1); %core edge noPref unknown nonabundant
        SpeciesEEresp(1).Category_prob = zeros(nbspecs,1); %highest posterior
        SpeciesEEresp(1).Category_prob_matrix = zeros(nbspecs,nbmaincats); %matrix of posterior probability (highest in each category)
        SpeciesEEresp(1).mean_abundance = zeros(nbspecs,1); %vector 1*nb species  % average of smoothed valid abundances 
        SpeciesEEresp(1).median_abundance = zeros(nbspecs,1);  %vector 1*nb species  % median of smoothed valid abundances 
        SpeciesEEresp(1).DEI = 0;  %DEI value
        SpeciesEEresp(1).DEIoptind = indDEIopt;  %index k of optimum DEI in [SpeciesEEresp(1).DEI … SpeciesEEresp(k).DEI  …SpeciesEEresp(end).DEI ]
        SpeciesEEresp(1).nbCategories_acrossDEIs = zeros(nbspecs,1); %1 value per species, same for all DEIs. Nb of categories given for all DEIs. Between 1 and 4.
        SpeciesEEresp(1).Fragmentation_Impact = zeros(nbspecs,1); %total abundance in landscape pixels / total abundance without EE in landscape pixels (from ab at each PC)
        SpeciesEEresp(1).Fragmentation_Impact_nowater = zeros(nbspecs,1); %total abundance in landscape pixels no water / total abundance without EE in landscape pixels no water (from ab at each PC)
        SpeciesEEresp(1).is_indicator = zeros(nbspecs,1); %1 value per species, computed for DEIopt only. boolean
        SpeciesEEresp(1).is_usableFIana = zeros(nbspecs,1); %1 value per species, computed for DEIopt only. boolean. if species has same category for DEIopt and DEI 1km
        SpeciesEEresp(1).EIsensitivity = zeros(nbspecs,1); %1 value per species, abundance on graph / ab on graph if species was insensitive to edge effects
        SpeciesEEresp(1).ind_AbPatRed_closestpattern = zeros(nbspecs,1); %indice in AbPatRef of the closest pattern to the species abundance curve
        SpeciesEEresp(2:length(LC))  = SpeciesEEresp(1);

               

        % ----------------------------------------------------------------------------------loop on DEI value
        for indDEI = 1 : nbDEIs %loop on DEI values

            %list of ScaleCategory names of potential patterns for this dataset and DEI
            potentialpattern_ScaleCategory = {AbPatRef(pattern_ispotential(:,indDEI)).ScaleCategory}'; %#ok<NODEF>
            %list of Category names of potential patterns for this dataset and DEI
            potentialpattern_Category = {AbPatRef(pattern_ispotential(:,indDEI)).Category}';
            %unique potential ScaleCategories and indices
            [~, ui] = unique(potentialpattern_ScaleCategory,'stable'); %list of potential patterns
            %Category names corresponding to the unique ScaleCategories (cannot find those with regex (eg Forest Ds Core .. and Forest Core)
            uniq_potentialpattern_ScaleCategory_corresCategory = potentialpattern_Category(ui);
           
            %HDmap
            HDmap = HD(indDEI).map;

            for indspec = 1:nbspecs %-------------------------------------------------------loop on species


                %extract valid measurements position:
                valid_mes = RespPCHD(indDEI).valid_mes_matrix(:,indspec);

                %smooth abundance PC/HD, valid measurements, for species indspec and DEI indDEI           
                speAb_smoo_valid = RespPCHD(indDEI).spemat_smoothPCHD(valid_mes,indspec);
                speAb_smoo_valid = speAb_smoo_valid .* double(speAb_smoo_valid > 0); % set negative values to zero

                %extract PC and HD of each plot (this has to be done per species because some plots are outliers and that depends on species)
                validplotPC = RespPCHD(indDEI).plot_PointCover(valid_mes);
                validplotHD = RespPCHD(indDEI).plot_HabitatDisturbance(valid_mes);
                
                
                %check if species has some abundance - if not no need to process or create figure
                if any(speAb_smoo_valid) %if any non zero abundance values, process

                    %--------------------------------- 1- INTERPOLATION, 'natural'

                    %if all plots have same PC (PID85), colinear, add -1 0 1
                    if length(unique(validplotPC))==1
                        validplotPC = validplotPC + round(rand(length(validplotPC),1)*2)-1;
                    end
                    if length(unique(round(validplotHD)))<=2 %if all the local covers are similar (PID 126 -> 130)
                        validplotHD = validplotHD + round(rand(length(validplotHD),1)*2)-1;
                    end

                    warning('off');
                    F = TriScatteredInterp(validplotPC,validplotHD,speAb_smoo_valid,'natural'); %#ok<*DTRIINT> %2012 
                    %F = scatteredInterpolant(validplotPC,validplotHD,speAb_smoo_valid,'natural','nearest');  %2015 
                    [Xfull, Yfull] = meshgrid(0:100,-100:100); %all combinations of PC (0 to 100) and EI (-100 to 100), integers
                    interpAbPCHD = F(Xfull(:), Yfull(:)); %get interpolated abondances on the full PC HD graph

                    %--------------------------------- 2- EXTRAPOLATION, 'nearest' (safest option)
                    posinterpvals = ~isnan(interpAbPCHD(:)); %the interpolated values are the non NaNs in interpAbPCHD, extrapolate around
                    F = TriScatteredInterp(Xfull(posinterpvals),Yfull(posinterpvals),interpAbPCHD(posinterpvals),'nearest'); 
                    %Use nearest-neighbor interpolation to replace these NaN values with the value at the nearest sample point
                    F.Method = 'nearest';
                    extrapAbPCHD = interpAbPCHD; 
                    extrapAbPCHD(~posinterpvals) = F(Xfull(~posinterpvals), Yfull(~posinterpvals));  % Use linear indexing to find the indices of the query points that produce NaN values when interpolated
                    extrapAbPCHD = extrapAbPCHD .* double(extrapAbPCHD > 0); %extrapAbPCHD is a full vector with no NaNs, and no neg values

                    %reshape on 201x101 graph
                    [graph_nrow, graph_ncol] = size(AbPatRef(1).matrix);
                    extrapAbPCHD=reshape(extrapAbPCHD,[graph_nrow, graph_ncol]);

                    %--------------------------------- 3- find the X points around each valid plots
                    plotmask = zeros(graph_nrow, graph_ncol);
                    indPlotHD = round(validplotHD)+101;
                    indPlotPC = round(validplotPC)+1;
                    indPlot_ongraph = sub2ind([graph_nrow, graph_ncol], indPlotHD , indPlotPC); %line indices
                    plotmask(indPlot_ongraph) = 1; % binary image where valid plot location is white

                    SE = strel('disk',2);
                    extrapzone = imdilate(plotmask,SE); %dilate around valid plots, 5*5, 25 points
                    extrapzone = logical(extrapzone); %make it boolean
                    
                    %--------------------------------- 4- PC, HD and abundance in extrapolated zone
                    extrap_speAb_smoo_valid = extrapAbPCHD(extrapzone(:)); %as a vector


                    %--------------------------------- 5- normalise abundance values
                    max_extrap_speAb_smoo_valid = max(extrap_speAb_smoo_valid); 
                    if max_extrap_speAb_smoo_valid > 0
                        extrap_speAb_smoonorm_valid = extrap_speAb_smoo_valid / max_extrap_speAb_smoo_valid;
                    else %avoid dividing by zero
                        extrap_speAb_smoonorm_valid = extrap_speAb_smoo_valid;
                    end

                    %--------------------------------- 6- compile training set : values of potential patterns for all points in extrap zone
                    indpotPat = 1; 
                    nbextrappoints = sum(extrapzone(:));
                    patternAb_extrap = zeros(nbextrappoints,sum(pattern_ispotential(:,indDEI)));
                    for indPat = 1:nbpatterns
                        if  pattern_ispotential(indPat,indDEI)
                            patternAb_extrap(:,indpotPat) = AbPatRef(indPat).matrix(extrapzone(:));
                            indpotPat = indpotPat + 1;
                        end
                    end

                    %--------------------------------- 7- Naive Bayes classification
                    %with smoothed data points
                    %[predicted_category, ~, posterior_vec] = classify(speAb_smoonorm_valid',patternAb_validplots_potential',potentialpattern_category,'diaglinear');    
                    %with values interpolated on PCLC graph :
                    [predicted_category, ~, posterior_vec] = classify( ...
                        extrap_speAb_smoonorm_valid', ...
                        patternAb_extrap', ...
                        potentialpattern_ScaleCategory,'diaglinear'); 


                    %--------------------------------- 8- Store category and category properties 

                    ispat_inpredcat = arrayfun(@(x) strcmp(x.ScaleCategory, predicted_category) , AbPatRef); %find indices of selected patterns in AbPatRef
                    ind_predicted_category = find(ispat_inpredcat);
                    ind_predicted_category = ind_predicted_category(1); %just take the first one, they are all in same category, we just extract names
            
                    %fill in output structure SpeciesEEresp:
                    SpeciesEEresp(indDEI).Category{indspec} = AbPatRef(ind_predicted_category).Category; %short name of category
                    SpeciesEEresp(indDEI).Habitat{indspec} = AbPatRef(ind_predicted_category).Habitat;
                    SpeciesEEresp(indDEI).EIPref{indspec} = AbPatRef(ind_predicted_category).EIPref;
                    SpeciesEEresp(indDEI).Category_prob(indspec) = max(posterior_vec); %highest posterior / just an indicator
                    SpeciesEEresp(indDEI).mean_abundance(indspec) = mean(speAb_smoo_valid); %vector 1*nb species  % average of smoothed valid abundances 
                    SpeciesEEresp(indDEI).median_abundance(indspec) = median(speAb_smoo_valid);  %vector 1*nb species  % median of smoothed valid abundances 
                    SpeciesEEresp(indDEI).DEI = LC(indDEI).DEI;  %DEI value                    
                    
                    %--------------------------------- find and store closest pattern to speab
                    %find closest pattern to current speab, amongst valid patterns AND patterns in attributed category
                    indpat_potential_inpredcat = find(pattern_ispotential(:,indDEI) .* ispat_inpredcat);
                    %compute sum of squared difference (on extrap zone)
                    sumsqd = arrayfun(@(x) sum((x.matrix(extrapzone(:)) - extrap_speAb_smoonorm_valid).^2), AbPatRef(indpat_potential_inpredcat));
                    %find indice of min simsqd
                    [~, indmin] = min(sumsqd);
                    %find indice of closest pattern in AbPatRef
                    SpeciesEEresp(indDEI).ind_AbPatRed_closestpattern(indspec) = indpat_potential_inpredcat(indmin);

                    %--------------------------------- 9- Store posterior probabilities for all main categories

                    %for each response category (forest core, forest edge, forest noPref, etc..) :
                    for indmaincat = 1:nbmaincats
                        
                        %indices of corresponding potential categories
                        boolcorrescats = strcmp(main_category_list{indmaincat}, uniq_potentialpattern_ScaleCategory_corresCategory);

                        %all posterior probablities in this main category
                        all_post_prob_incat = posterior_vec(boolcorrescats);

                        if isempty(all_post_prob_incat) %if category could not be assessed (no potential paterns in there)
                            SpeciesEEresp(indDEI).Category_prob_matrix(indspec,indmaincat) = 0; %posterior prob is zero
                        else                
                            SpeciesEEresp(indDEI).Category_prob_matrix(indspec,indmaincat) = sum(all_post_prob_incat); %posterior prob of main category as sum of all sub categories
                        end
                    end


                    %--------------------------------- 10- Compute fragmentation impact
                    %EIPref = SpeciesEEresp(indDEI).EIPref{indspec};
                    %if ~strcmp(EIPref,'Unknown') && ~strcmp(EIPref,'Nonabundant') && ~strcmp(EIPref,'Nonsensitive')
                    %if the species is sensitive to edge effects -----> do anyway.

                    %---- not using full extrapolation - but assume 0 value where people have not measured

                    %redilate plotmask, but much more
                    SE2 = strel('disk',8);
                    extrapzone2 = imdilate(plotmask,SE2); %dilate around valid plots
                    
                    %extrapolated abundance on graph around measured plots dilated 8:
                    extrapAbPCHD_onextrapzone = extrapAbPCHD .* extrapzone2; %extrap val in zone, 0 outside
           

                    %----transposing abundance values on landscape

                    %(line) indices of the PC/HD pairs from map on the 201*101 pchd graph  | HD is in rows
                    valpixHD = round(HDmap(:))+101;
                    valpixPC = round(PCmap(:))+1;
                    indpix_ongraph = sub2ind([graph_nrow, graph_ncol], valpixHD , valpixPC); %line indices

                    %extrap abundance values for all map pixels
                    extrapAb_onmap = extrapAbPCHD_onextrapzone(indpix_ongraph); %using zero mod - full extrapolation
                    extrapAb_onmap = reshape(extrapAb_onmap, size(PCmap)); %reshape same size as map
                    %extrapAb_onmap contains the abundance of the species in map format, for all measured value of PC and EI
     
                    extrapzone_onmap = extrapzone2(indpix_ongraph); %map zone for which we can extrapolate abundance values (1=extrap ok, 0=unknown)
               
                    %----computing abundance on map without edge effects

                    %max of each PCHC graph column (i.e. max abundance for each PC)
                    maxAbforeachPC = max(extrapAbPCHD_onextrapzone); %using zero mod full extrapolation
                    %--- create nonsensi graph, then port on map
                    maxAbforeachPC = repmat(maxAbforeachPC,graph_nrow,1);
                    maxAbforeachPC_onextrapzone = maxAbforeachPC .* extrapzone2; %graph of ab values, no variations along EI axis, in extrap zone

                    %abundance extrapolated on measured zone on map, from PC sensitivity
                    nonsensi_extrapAb_onmap = maxAbforeachPC_onextrapzone(indpix_ongraph); 
                    nonsensi_extrapAb_onmap = reshape(nonsensi_extrapAb_onmap, size(PCmap)); %reshape same size as map
                  
                    %---- compute 1 - ratio :
                    Fragmentation_Impact = 1 - sum(extrapAb_onmap(:)) / sum(nonsensi_extrapAb_onmap(:));

                    %---- compute 1 - ratio no water :
                    %counting 0 abundance where there is water
                    Fragmentation_Impact_nowater = 1 - sum(extrapAb_onmap(:) .* ~mask.map(:)) / sum(nonsensi_extrapAb_onmap(:) .* ~mask.map(:));
                    %impact of edge sensitivity, if ratio close to 1 , no impact.
                    
                    %Fragmentation impact is the ratio of the sum of abundance on measured zone on map, 
                    %divided by the sum of abundance nonsensi to EI on measured on map. 
                    %Excluding map zones where TC/EI combinations have not been measured.
                    

                    %---- store in SpeciesEEresp
                    SpeciesEEresp(indDEI).Fragmentation_Impact(indspec) = Fragmentation_Impact; %#ok<*AGROW>
                    SpeciesEEresp(indDEI).Fragmentation_Impact_nowater(indspec) = Fragmentation_Impact_nowater;
                        
                    
                    %--------------------------------- 11- Compute Edge Influence Sensitivity
                    
                    %EI sensitivity is the ratio of abundance (EIsensi and not EI sensitive) on the graph, in order to compare species measured in different landscapes
                    
                    %---- computing the impossible zone on graph: (EI range is 100 for all TC), this is not extrap zone
                    [tcvals, eivals] = meshgrid(0:100,-100:100);
                    validzone = ones(201,101);
                    validzone(eivals>100-tcvals)=0; %EI cannot be higher than 100-TC
                    validzone(eivals<-tcvals)=0; %EI cannot be lower than -TC
                    
                    %---- full graph extrapolation of abundance from TC and EI  
                    abEITC = extrapAbPCHD .* validzone; %extrapolated abundance on full graph, force 0 in invalid zone
                    %---- full graph extrapolation o0f abundance from TC only
                    maxAbforeachPC_fullgraph = repmat(max(extrapAbPCHD .* validzone),graph_nrow,1); %grab the max of each column in valid zone of full graph of predicted abundance, then replicate the row
                    abTC = maxAbforeachPC_fullgraph .* validzone; %max of each PCHC full graph validzone column (i.e. max abundance for each PC), force 0 in invalid zone  
                    %---- compute sensitivity, ratio of the two graphs (measured and PC estimated:
                    EIsensitivity = 1 - sum(abEITC(:))/sum(abTC(:)); %ratio of ab on full graph  ; and ab insensitive to EI on full graph ==> EI range
                    SpeciesEEresp(indDEI).EIsensitivity(indspec) = EIsensitivity; %---- store in SpeciesEEresp
                    
                    
                    %end % end on edge sensitivity condition 




                    %--------------------------------------------------------------------------------------------------------------------------
                    %----------------- MAKE THE FIGURE ----------------------------------------------------------------------------------------

                    %for each species save a figure showing:
                    %1/species abundance in valid plot on PCmap + give species number, category, name, EIsensi, and FI
                    %2/species abundance in valid plot on EImap + give DEI used
                    %3/species abundance in valid plot on PCEI graph + give species subcategory
                    %4/species abundance extrapolated on graph
                    %5/species abundance extrapolated on graph from PC only
                    %6/species abundance in valid plot in Forest only and Matrix only
                    %7/species abundance extrapolated on map with EI and PC
                    %8/species abundance extrapolated on map with PC only (if species was insensitive)
                    %9/species abundance in valid plots wrt Distance

                    h=figure('visible','off','color','w','position',[20 20 1280 1024]);
                    cop = [0.7 0.7 0.7; jet(101)];
                    fss = 6;
                    fsh = 7;
                    
                    %compute a 'no color'/outside of extrap zone value, 1% of max ab
                    maxAb = max(speAb_smoo_valid(:));
                    NANval = -maxAb/100;

                    
                    %also extract plot coordinates on map
                    plot_row = [plotLoc.pR]'; plot_row = plot_row(valid_mes);
                    plot_col = [plotLoc.pC]'; plot_col = plot_col(valid_mes);

                    %and distance to nearest edge + smooth abundances
                    plot_dist = [plotLoc.dist]'; plot_dist = plot_dist(valid_mes);
                    abdist = species.matrix(valid_mes,indspec); %raw species abundance
                    [plot_dist, od] = sort(plot_dist); abdist = abdist(od);
                    propplottosmooth = 0.25;
                    abdist_smoo = smooth(plot_dist, abdist, propplottosmooth, 'lowess');


                    %1/----------------------------------------------species abundance in valid plot on PCmap
                    subplot(2,3,1);
                    PCmap_disp = PC.map; PCmap_disp(1,1)=-1; %to avoid displaying 0 in gray
                    imagesc(PCmap_disp); colormap(cop); axis image; ch = colorbar('southoutside'); 
                    hold on
                    scatter(plot_col, ...
                    plot_row, ...                
                    sqrt(speAb_smoo_valid)*100+1, ...
                    'k*');  
                    set(gca,'Xtick',[]);set(gca,'Ytick',[]);
                    set(gca,'fontsize',fss);
set(get(ch,'XLabel'),'String','Point Cover (PC)','fontsize',fss); colorbar_pos = get(ch,'Position'); axes_pos = get(gca, 'Position');
new_colorbar_pos = colorbar_pos; new_colorbar_pos(4) = 0.5*colorbar_pos(4); new_colorbar_pos(2) = axes_pos(2)-axes_pos(4)/5;
set(ch,'Position',new_colorbar_pos); set(gca,'Position',axes_pos); %this move colorbar depending on axes position and height
                    title({['#' num2str(indspec) ': ' species.names{indspec} ' => ' SpeciesEEresp(indDEI).Category{indspec}], ...
                        ['EI sensitivity = ' num2str(round(SpeciesEEresp(indDEI).EIsensitivity(indspec)*100)) ' %' ...
                        ', FI = ' num2str(round(SpeciesEEresp(indDEI).Fragmentation_Impact_nowater(indspec)*100)) ' %'], ...
                        'Abundance on PC map'}, ...
                    'fontsize',fsh);                    

                    %2/----------------------------------------------species abundance in valid plot on EImap
                    subplot(2,3,2);
                    HDmap_disp = HD(indDEI).map; HDmap_disp(1,1)=-101; %to avoid displaying -100 in gray
                    imagesc(HDmap_disp); colormap(cop); axis image;  ch = colorbar('southoutside'); 
                    hold on
                    scatter(plot_col, ...
                    plot_row, ...                
                    sqrt(speAb_smoo_valid)*100+1, ...
                    'k*');  
                    set(gca,'Xtick',[]);set(gca,'Ytick',[]);
                    set(gca,'fontsize',fss);
set(get(ch,'XLabel'),'String','Edge Influence (EI)','fontsize',fss); colorbar_pos = get(ch,'Position'); axes_pos = get(gca, 'Position');
new_colorbar_pos = colorbar_pos; new_colorbar_pos(4) = 0.5*colorbar_pos(4); new_colorbar_pos(2) = axes_pos(2)-axes_pos(4)/5;
set(ch,'Position',new_colorbar_pos); set(gca,'Position',axes_pos);
                    title({['maxDEI=' num2str(LC(indDEI).DEI) 'm'], ...
                        ' ', ...
                        'Abundance on EI map'}, ...
                    'fontsize',fsh);
                
                    %3/----------------------------------------------species abundance in valid plot on PCEI graph
                    subplot(2,3,3);
                    scatter(validplotPC, ...
                    validplotHD, ...                
                    speAb_smoo_valid*50+1, ...
                    speAb_smoo_valid);  colormap(cop); colorbar;
                    set(gca,'fontsize',fss);
                    ylabel('Edge Influence (EI)'); xlabel('Point Cover (PC)');  
                    title({['=> ' SpeciesEEresp(indDEI).Habitat{indspec} ' ' SpeciesEEresp(indDEI).EIPref{indspec}], ...
                        'Abundance wrt PC and EI'}, ...
                    'fontsize',fsh);                
                    hold on; plot([0 100],[0 0],'k'); plot([0 100],[50 0],'k--'); plot([0 100],[0 -50],'k--'); 
                    hold off; axis image; axis([0 100 -100 100]);
                                      
                    %4/----------------------------------------------species abundance in valid plot in Forest only and Matrix only
                    subplot(2,3,4);
                    indForest = find(validplotPC > 65);
                    xf = abs(validplotHD(indForest)); yf = speAb_smoo_valid(indForest); ysf=smooth(xf,yf,propplottosmooth,'lowess');
                    [xf, nx] = sort(xf); yf=yf(nx); ysf=ysf(nx);
                    plot(xf, ...
                    yf, ...
                    'xg','linewidth',2);
                    hold on;  
                    indMatrix = find(validplotPC < 35);
                    xm = abs(validplotHD(indMatrix)); ym = speAb_smoo_valid(indMatrix); ysm=smooth(xm,ym,propplottosmooth,'lowess');
                    [xm, nx] = sort(xm); ym=ym(nx); ysm=ysm(nx);
                    plot(xm, ...
                    ym, ...
                    'xr','linewidth',2);  
                    plot(xf,ysf,'-g','linewidth',2);
                    plot(xm,ysm,'-r','linewidth',2);
                    set(gca,'fontsize',fss);
                    xlabel('| Edge Influence (EI) |'); ylabel('Species abundance');  
                    title({'Abundance vs |EI|'}, ...
                    'fontsize',fsh);                
                    plot([50 50],[0 maxAb],'k--');hold off
                    legend('in forest','in matrix');
                    
                    %5/----------------------------------------------species abundance in valid plots wrt Distance
                    subplot(2,3,5);
                    plot(plot_dist,abdist,'ob');
                    hold on
                    plot(plot_dist,abdist_smoo,'rx-');
                    set(gca,'fontsize',fss);
                    ylabel('Species abundance'); xlabel('Distance to nearest edge');  
                    title('Abundance vs Distance to nearest edge', ...
                    'fontsize',fsh);                
                    hold off; 
                
             
                    %6/----------------------------------------------species abundance extrapolated on map with EI and PC
                    subplot(2,3,6);
                    tp = extrapAb_onmap.* ~mask.map;
                    tp(extrapzone_onmap == 0) = NANval;
                    imagesc(tp); colormap(cop); axis image;  ch = colorbar('southoutside');
                    set(gca,'Xtick',[]);set(gca,'Ytick',[]);
                    set(gca,'fontsize',fss);
set(get(ch,'XLabel'),'String','predicted abundance','fontsize',fss); colorbar_pos = get(ch,'Position'); axes_pos = get(gca, 'Position');
new_colorbar_pos = colorbar_pos; new_colorbar_pos(4) = 0.5*colorbar_pos(4); new_colorbar_pos(2) = axes_pos(2)-axes_pos(4)/5;
set(ch,'Position',new_colorbar_pos); set(gca,'Position',axes_pos);
                    title({'Abundance extrapolated from PC and EI'}, ...
                        'fontsize',fsh);


                    %                     %8/----------------------------------------------species abundance extrapolated on map with PC only (if species was insensitive)
                    %                     subplot(3,3,8);
                    %                     tp = nonsensi_extrapAb_onmap.* ~mask.map;
                    %                     tp(extrapzone_onmap == 0) = NANval;
                    %                     imagesc(tp); colormap(cop); axis image;  ch = colorbar('southoutside'); 
                    %                     set(gca,'Xtick',[]);set(gca,'Ytick',[]);
                    %                     set(gca,'fontsize',fss);
                    % set(get(ch,'XLabel'),'String','predicted abundance','fontsize',fss); colorbar_pos = get(ch,'Position'); axes_pos = get(gca, 'Position');
                    % new_colorbar_pos = colorbar_pos; new_colorbar_pos(4) = 0.5*colorbar_pos(4); new_colorbar_pos(2) = axes_pos(2)-axes_pos(4)/1.9;
                    % set(ch,'Position',new_colorbar_pos); set(gca,'Position',axes_pos);
                    %                     title({'Extrapolated Abundance wrt PC only'}, ...
                    %                     'fontsize',fsh); 
                
                

                    saveas(h,fullfile(runvar.folder_Data_Output,[runvar.species_name '_' num2str(indspec) '_DEI_' num2str(runvar.DEI) '.tif']));
                    close(h);
                    
                    
                else %if species is absent
                    %fill in output structure SpeciesEEresp:
                    SpeciesEEresp(indDEI).Category{indspec} = 'Absent'; %short name of category
                    SpeciesEEresp(indDEI).Habitat{indspec} = 'Absent';
                    SpeciesEEresp(indDEI).EIPref{indspec} = 'Absent';
                    SpeciesEEresp(indDEI).Category_prob(indspec) = NaN; %highest posterior
                    SpeciesEEresp(indDEI).mean_abundance(indspec) = 0; %vector 1*nb species  % average of smoothed valid abundances 
                    SpeciesEEresp(indDEI).median_abundance(indspec) = 0;  %vector 1*nb species  % median of smoothed valid abundances 
                    SpeciesEEresp(indDEI).DEI = LC(indDEI).DEI;  %DEI value
                    SpeciesEEresp(indDEI).Fragmentation_Impact(indspec) = NaN; %#ok<*AGROW>
                    SpeciesEEresp(indDEI).Fragmentation_Impact_nowater(indspec) = NaN;
                    SpeciesEEresp(indDEI).EIsensitivity(indspec) = NaN;
                    
                    %no fig
                    
                end %end of condition on nonzero abundance values

            end %end loop on species
        end %end loop on DEI values

       




        %-------------------store .mats
        save(runvar.species_category_matfile ,'SpeciesEEresp');

        msg=['computed species category and FI for ' runvar.species_name]; dispwrite_log(runvar, param, msg)

    else
        msg=['previously computed: species category and FI ' runvar.species_category_matfile]; dispwrite_log(runvar, param, msg)

    end %end of condition on output existence

end %end function

    

            
            
            
            
            
  
