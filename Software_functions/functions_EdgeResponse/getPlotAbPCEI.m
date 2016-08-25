% getPlotAbPCEI
% get Point cover and Habitat Disturbance values for each plot and save in plotLoc, for opt DEI (just for 1 DEI)
% for each species remove abundance outliers 
% smooth abundance values wrt PC and HD and store as mat file (with valid positions)

function runvar = getPlotAbPCEI(runvar, param)


    %INPUT VALUE--------------------------------------------------------------

    %nb of plot to smooth in the PCHD graph (span / window)
    %nbplotstosmooth = 51;
    %proportion of plots to smooth in the PCHD graph (span):
    propplottosmooth = 0.25;  %best to use a fixed prorportion as opposed to a number of plots. because interval is bonded. because huge variations in number of plots across datasets.

    %create filename for RespPCHD structure
    runvar.RespPCHD_matfile = fullfile(runvar.folder_Data_Matlab, [runvar.species_name '_DEI_' num2str(runvar.DEI) param.RespPCHD]);
    
    %-------------------------------compute PC and EI plots values, and smooth species abundances at plot
    if ~exist(runvar.RespPCHD_matfile,'file')
    
        %load PCmap, plotLoc and EI map
        load(runvar.PCmap_matfile);
        load(runvar.plot_matfile);
        load(runvar.LCmap_matfile);
    
        %------- extract Point Cover and plot distances, and all DEI used

        %Point Cover (PointCover):    
        PointCover = PC.map;
        
        %DEI used
        indDEIopt = 1;
        DEIopt = runvar.DEI;


        %------------------------- Part 1: update plotLoc with Point Cover and Local Cover, HD, DEI opt for all plots

        %--------------------- extract PC and LC and HD for all plots and save 

        %number of plots
        numplots = length(plotLoc);

        %plot image indices: convert pR and pC to single index of plots
        pR = [plotLoc.pR]';
        pC = [plotLoc.pC]';
        indPlot_onmap = sub2ind(size(PointCover),pR,pC);  %assuming all plots do fall within image (no NaN)
        mapcoord_match = ~any(isnan(indPlot_onmap(:)));
        
        %if not all plots fall on map provided send error:
        %prepare to throw error if needed:
        var_to_assert = mapcoord_match;
        erreason = ['not all census points coordinate fall within the extent of the map provided: ' runvar.PCmap_name];
        errID = 'inputs';
        check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if plots fall outside the map
        
        % plot Point Cover:
        plot_PointCover = PointCover(indPlot_onmap);

        %extracting Local Cover for optimal DEI:
        LocalCover_opt = LC(indDEIopt).map;

        %extracting habitat disturbance for optimal DEI:
        HD_opt = HD(indDEIopt).map; 

        % plot Local Cover for optimum DEI:
        plot_LCopt = LocalCover_opt(indPlot_onmap);
        % plot Habitat Disturbance for optimum DEI:
        plot_HDopt = HD_opt(indPlot_onmap);


        if ~isfield(plotLoc,'HabitatDisturbanceOpt') %if not already entered (check needed for presence/absence datasets

            %---------------save plotPC and plotLCopt, plotHDopt, in plotLoc

            %add plot pointcover in plotLoc structure:

            cellarr_plotPC = num2cell(plot_PointCover(:));

            [plotLoc(1:numplots).PointCover] = cellarr_plotPC{:};

            %add plot Local cover for optimal DEI in plotLoc structure:

            cellarr_plotLC = num2cell(plot_LCopt(:));

            [plotLoc(1:numplots).LocalCoverOpt] = cellarr_plotLC{:};

            %add plot Habitat Disturbance for optimal DEI in plotLoc structure:

            cellarr_plotHD = num2cell(plot_HDopt(:));

            [plotLoc(1:numplots).HabitatDisturbanceOpt] = cellarr_plotHD{:};

            %update plotLoc file
            save(runvar.plot_matfile,'plotLoc');

            msg=['plot PC and LCopt and HDopt added in plotLoc for ' runvar.plot_matfile]; dispwrite_log(runvar, param, msg)
        end



        %------------------------- Part 2: compute smooth abundance w.r.t. PC and HD (for all DEI, for all species)

        msg = ['... smoothing each species abundance wrt PC/EI for ' runvar.species_name]; dispwrite_log(runvar, param, msg)
        
        % suitable for abundance or normalised abundance datasets

        %load species file
        load(runvar.species_matfile);


        %extract species matrix and species names
        spemat = species.matrix; %spcies matrix (rows ordered same as plotLoc)
        nbspec = size(spemat,2); %number of species (columns)

        %number of DEI for which LC was computed
        nbDEI = length(LC);

        %RespPCHD structure - memory allocation
        RespPCHD(1).spemat_smoothPCHD = NaN(numplots,nbspec);
        RespPCHD(1).valid_mes_matrix = zeros(numplots,nbspec);
        RespPCHD(1).plot_PointCover = zeros(numplots,1);
        RespPCHD(1).plot_LocalCover = zeros(numplots,1);
        RespPCHD(1).plot_HabitatDisturbance = zeros(numplots,1);
        RespPCHD(1).DEI = 0;
        RespPCHD(1).DEIopt = 0;
        RespPCHD(1).indDEIopt = 0;
        RespPCHD(2:nbDEI) =  RespPCHD(1); %copy over


        for indDEI = 1:nbDEI  %for each DEI 

            %--------------------------- PC/HD at plot locations

            %extract Local Cover computed with LC(indDEI).DEI at plot locations
            plot_LocalCover = LC(indDEI).map(indPlot_onmap);
            %extract Habitat Disturbance computed with LC(indDEI).DEI at plot locations
            plot_HabitatDisturbance = HD(indDEI).map(indPlot_onmap);

            %--------------------------- remove abundance outliers and smooth over PC|HD, species by species

            %smooth abundance matrix - allocation
            spemat_smoothPCHD = NaN(numplots,nbspec);
            %valid measurement matrix - allocation
            valid_mes_matrix = false(numplots,nbspec);


            %for all species in spemat
            for indspec=1:nbspec   

                %abundance vector for species indspec:
                abvals = spemat(:,indspec);

                %initialize valid measurement to true (numplots measurements)
                valid_mes = true(numplots,1); %valid measurements


                %--- removing outliers (plots where the difference of their abundance with mean abundance is larger than 10 times the std of abundance)
                multconst = 10;
                indOutliers = abs(abvals - mean(abvals)) > multconst * std(abvals); %is equal to 1 when plot is outlier

                if any(indOutliers)  %if there are some outliers, set their value to false in valid_mes
                    valid_mes = ~indOutliers; %is equal to false where plot is outlier
                    %disp(['warning species #' num2str(indspec) ', ' spnames{indspec} ': ' num2str(sum(indOutliers)) ' outlier(s) removed']);
                end


                %--- smooth abundance values in PC / HD graph
                x = plot_PointCover(valid_mes);         % x plot point cover without outliers
                y = plot_HabitatDisturbance(valid_mes); % y plot habitat disturbance without outliers
                z = abvals(valid_mes);                  % z plot abundance without outliers

                f=fit([x,y],z, 'lowess', 'Span', propplottosmooth);       % build interpolant smooth - smooth 25% points
                abvals_PC_HD = feval(f,[x,y]);  % evaluate interpolant at plot locations on graph PCHD

                if any(isnan(abvals_PC_HD)) %if some NaNs it means that more than nb of smoothed plots have the same PCHD
                    f=fit([x,y],z, 'lowess', 'Span', propplottosmooth * 2);       % smooth by 50%, that should fix it
                    abvals_PC_HD = feval(f,[x,y]);  % evaluate interpolant at plot locations on graph PCHD
                    disp('warning : more than 25% of plots with same PC and HD'); 
                end

                abvals_PC_HD = abvals_PC_HD.*double(abvals_PC_HD>0); % set negative values to zero


                %--- store smooth abundance and valid measurements in matrices
                spemat_smoothPCHD(valid_mes,indspec) = abvals_PC_HD; %plot with invalid measurements are NaN
                valid_mes_matrix(:,indspec) = valid_mes;

            end %end loop on species indspec

            %store smooth abundance matrix, valid_mes matrix, plot point cover, LC and HD values at plot locations, and DEI value in structure
            RespPCHD(indDEI).spemat_smoothPCHD = spemat_smoothPCHD;
            RespPCHD(indDEI).valid_mes_matrix = valid_mes_matrix;
            RespPCHD(indDEI).plot_PointCover = plot_PointCover;     %store it unrounded
            RespPCHD(indDEI).plot_LocalCover = plot_LocalCover;     %store it unrounded
            RespPCHD(indDEI).plot_HabitatDisturbance = plot_HabitatDisturbance;     %store it unrounded
            RespPCHD(indDEI).DEI = LC(indDEI).DEI;
            RespPCHD(indDEI).DEIopt = DEIopt;
            RespPCHD(indDEI).indDEIopt = indDEIopt;

        end %end loop on DEI indDEI

        %save RespPCHD as mat file for PID i
        save(runvar.RespPCHD_matfile,'RespPCHD');

        msg=['species abundance smoothed wrt PC/EI for ' runvar.species_name]; dispwrite_log(runvar, param, msg)
        
    else

        msg=['previously computed: species smooth abundance ' runvar.RespPCHD_matfile]; dispwrite_log(runvar, param, msg)

    end %end of condition on output existence




end
    
    
    
     
    
    
    
    