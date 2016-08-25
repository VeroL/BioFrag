% function projectPlotsLatLon

% project the plot from Lat - Lon to x - y and indices of rows and columns
% using projection of associated raster, and stored in ginfo and R

% stores projected coords and indices as new fields in the plotLoc structure 

% 1 x number of plots struct array with fields:
%     name
%     lat
%     lon
%     X
%     Y
%     pR
%     pC

function [plotLoc, mapcoord_match] = projectPlotsLatLon(plotLoc, PointCover)

    %number of plots
    numplots = length(plotLoc);
    
    %if fields X and Y are empty, project from lat|lon
    if isempty(plotLoc(1).X) %assuming that if X is empty then lat and lon are filled
        lat = [plotLoc.lat]';
        lon = [plotLoc.lon]';
    
        %project lat and lon of plot locations using given projection:
        [x, y] = projfwd(PointCover.ginfo,lat,lon);
        
        %add in plotLoc structure:
        
        %store projected coordinate into cell array, make sure each variable is in column
        projcoord_cellarr = num2cell([x(:) y(:)]); 
    
        [plotLoc(1:numplots).X] = projcoord_cellarr{:,1}; % projected map coord x
        [plotLoc(1:numplots).Y] = projcoord_cellarr{:,2}; % projected map coord y
    else
        %grab x and y coordinates from plotLoc
        
        x = [plotLoc.X]';
        y = [plotLoc.Y]';
    end
        
    
    %get pixel indices of plot locations in landscape map corresponding to input geo location (extent)
    [pR, pC]=PointCover.R.worldToSub(x, y);

    %add in plotLoc structure:
 
    %store projected coordinate into cell array, make sure each variable is in column
    projcoord_cellarr = num2cell([pR(:) pC(:)]); 
   
    [plotLoc(1:numplots).pR] = projcoord_cellarr{:,1}; % row in matrix map
    [plotLoc(1:numplots).pC] = projcoord_cellarr{:,2}; % column in matrix mat
    
    
    %check that plot row and col indices fall within the extent of the map provided:
    indPlot_onmap = sub2ind(size(PointCover.map),pR,pC);  %linear index of plots on map. NaN if exceeds map dimanesion.
    mapcoord_match = ~any(isnan(indPlot_onmap(:))); %true if no nan (all plots on map), 0 otherwise
    

end



