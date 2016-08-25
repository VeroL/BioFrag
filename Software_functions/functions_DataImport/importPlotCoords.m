% function importPlotCoords import plots from CSV file and store in a matlab structure array

% reads PID####_Plots.csv of format:
% plot name | lat | lon
%  string   | num | num

%converts first column to string if necessary
%order by natural order

% The plots coordinates are stored in the structure PlotLoc of format:
% 1 x number of plots struct array with fields:
%     name
%     lat
%     lon

% plotLoc(1)
% ans = 
%     name: 'B1'
%      lat: -12.18993333
%      lon: 44.23013333
%      
% {plotLoc.name}
% ans = 
%   Columns 1 through 12
%     'B1'    'B2'    'B3'    'B4'    'D1'    'D2'    'D3'    'D4'    'D5'    'D6'    'H1'    'H2'  ...
%     
% [[plotLoc.lat]' [plotLoc.lon]']
% ans =
%               -12.18993333               44.23013333
%               -12.18688333               44.22553333
%                  -12.18125               44.22446667
%                   -12.1762               44.22073333
%               -12.22306667               44.43508333
%               ...
              
% numplots = length(plotLoc);


function [plotLoc, headerline_ok] = importPlotCoords(pathPlotCsv)
    
    headerline_ok = true;
    
    %read the csv file - store into a cell array (each comma separatedvalue in individual cell
    [~,~,cellcsv] = xlsread(pathPlotCsv);
    
    %the format of the csv file is assumed to be either:
    % plot name | lat | lon
    % string    | num | num
    
    %or
    % plot name | x | y
    % string    | num | num
    
    %grab headerline:
    headerline = cellcsv(1,:);
    
    %find lat/lon or x/y, ignore case, look for 'lat' and 'lon' so that latitude and longitude also accepted
           
    %search for 'lat' (case insensitive), compare first 3 characters in all strings of headerline:
    lat_col = find(strncmpi('lat',headerline,3)); %col # if lat found, empty otherwise
    %search for 'lon':
    lon_col = find(strncmpi('lon',headerline,3)); %col # if lon found, empty otherwise
    
    %search for 'x', compare full string (case insensitive):
    x_col = find(strcmpi('x',headerline)); %col # if x found, empty otherwise
    %search for 'y':
    y_col = find(strcmpi('y',headerline)); %col # if y found, empty otherwise
   
    
    %check if the first column contains string, if not convert
    plotnames=cellcsv(2:end,1);
    if iscellstr(plotnames) == 0  %if not cell array of string
        plotnames = cellfun(@num2str, plotnames, 'UniformOutput',false);
    end
    
    %store number of plots
    numplots=length(plotnames);
    
    %build empty structure:
    plotLoc.name = [];
    plotLoc.lat = [];
    plotLoc.lon = [];
    plotLoc.X = [];
    plotLoc.Y = [];
    %copy over number of plots times
    plotLoc(2:numplots,:)=plotLoc(1);
    

    %fill in plot coordinates structure:    
    [plotLoc(1:numplots).name] = plotnames{:};
    
    if ~isempty(lat_col) && ~isempty(lon_col)  %if both lat and lon have been found:
        [plotLoc(1:numplots).lat] = cellcsv{2:end,lat_col};
        [plotLoc(1:numplots).lon] = cellcsv{2:end,lon_col};
    elseif ~isempty(x_col) && ~isempty(y_col) %if both x and y have been found:        
        [plotLoc(1:numplots).X] = cellcsv{2:end,x_col};
        [plotLoc(1:numplots).Y] = cellcsv{2:end,y_col};
    else %if lat|lon or x|y have not been found
        headerline_ok = false; %send back error if cannot find string pairs
    end
    
    
    %check if plot names are in natural order: a1 a2 a3 a10 a11...
    [~, indp]=sort_nat(plotnames);
    
    %if not reorder:
    if all(indp == (1:length(indp))') == 0 
        plotLoc = plotLoc(indp);    
    end

end