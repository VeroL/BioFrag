%getPlotDistances
%compute plot distances

function runvar = getPlotDistances(runvar, param)

    %load PCmap, plotLoc and distance map
    load(runvar.PCmap_matfile);
    load(runvar.plot_matfile);
    load(runvar.Distmap_matfile);
  
    if ~isfield(plotLoc,'dist') %#ok<NODEF> %if dist info not stored in plotLoc
        %--------------------- extract dist for all plots

        %distance map:
        dist_map = D2nE.map;

        %resolution
        res = PC.res;

        %plot image indices: convert pR and pC to single index of plots
        pR = [plotLoc.pR]';
        pC = [plotLoc.pC]';
        indPlot_onmap = sub2ind(size(dist_map),pR,pC);  %assuming all plots do fall within image (no NaN)
        mapcoord_match = ~any(isnan(indPlot_onmap(:)));
                
        %if not all plots fall on map provided send error:
        %prepare to throw error if needed:
        var_to_assert = mapcoord_match;
        erreason = ['not all census points coordinate fall within the extent of the map provided: ' runvar.PCmap_name];
        errID = 'inputs';
        check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if plots fall outside the map
        

        %number of plots
        numplots = numel(indPlot_onmap);

        % plot distance to nearest edge:
        plot_dist = dist_map(indPlot_onmap);
        %compute distance from pixel center rather than outer edge %consequence: res/2 or -res/2 is minimal distance to edge
        plot_dist = plot_dist - sign(plot_dist)*res/2; %remove res/2 if dist positive (in forest), add if dist is negative (outside forest)

        %add plot distance in plotLocClip structure:

        cellarr_plotdist = num2cell(plot_dist(:));

        [plotLoc(1:numplots).dist] = cellarr_plotdist{:};

        %update plotLoc file
        save(runvar.plot_matfile,'plotLoc');

        msg=['plot distances added for ' runvar.plot_matfile]; dispwrite_log(runvar, param, msg)
    else
        msg=['previously computed: plot distances ' runvar.plot_matfile]; dispwrite_log(runvar, param, msg)
    end
   
    
end
    
    