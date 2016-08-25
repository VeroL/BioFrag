% function importSpeciesMatrix and store in a structure

%the format of the csv file is assumed to be:
% plot      | species name 1 | species name 2 | ..
% plot name | abundance      | abundance      | ..
% string    | num            | num

%converts first column and first line to string if necessary
%order by natural order
%compute abundant species

% The csv file info are stored as 1 element structure species:
% PID####_species.mat
% species.matrix                   matrix, double, plot in rows, species in columns 
% species.plotNames                cell array of str, 1 col
% species.names                    cell array of str, 1 col
% species.isAbundant               boolean vector, 1 for abundant species, 1 * nb of species

%abundant species = species with more than 3 individuals in at least one plot when abundance is meaured
%or present in at least 3 plots when presence is measured

function species = importSpeciesMatrix(pathSpeciesCsv)

    %read the csv file - store into a cell array (each comma separatedvalue in individual cell
    [~,~,cellcsv] = xlsread(pathSpeciesCsv);
    
    %the format of the csv file is assumed to be:
    % plot      | species name 1 | species name 2 | ..
    % plot name | abundance      | abundance      | ..
    % string    | num            | num
    
        
    species_matrix = cell2mat(cellcsv(2:end,2:end));  %grab the species matrix from row 2, column 2
    
    %sometimes this bugs because there is a space in an empty cell. find
    %faulty cell and remove space manually:
    %     toto=cellfun(@(x) class(x), cellcsv(2:end,2:end),'UniformOutput',false);
    %     toto2 = cellfun(@(x) strcmp(x,'double'),toto);
    %     [i, j] = find(toto2==0)
    
    %this is to remove the NaNs from the cell array, however might not be needed.
    %     species_cell = cellcsv(2:end,2:end); %grab the species matrix from row 2, column 2
    %     indNaN = cellfun(@(x) isnan(x), species_cell);
    %     species_cell(indNaN) = {0}; %replace NaNs by 0s in the cell array
    %     species_matrix = cell2mat(species_cell);
   
    %remove NaNs - replace by zeros
    species_matrix(isnan(species_matrix)) = 0;
    
    %check if the first column (plot names) contains string, if not convert
    plotnames=cellcsv(2:end,1);
    if iscellstr(plotnames) == 0  %if not cell array of string
        plotnames = cellfun(@num2str, plotnames, 'UniformOutput',false);
    end
    
   
    %check if the first line (species names) contains string, if not convert
    speciesnames = cellcsv(1,2:end)'; %transpose because num2str works on column vectors
    if iscellstr(speciesnames) == 0
        speciesnames = cellfun(@num2str, speciesnames, 'UniformOutput',false);
    end
    
        
    %check if plot names are in natural order: a1 a2 a3 a10 a11...
    [~, indp]=sort_nat(plotnames);
            
    %if not reorder:
    if all(indp == (1:length(indp))') == 0 
        plotnames = plotnames(indp);  %the plot names
        species_matrix = species_matrix(indp,:); %the species matrix
    end
    
    ab_thresh = 2;
    if max(species_matrix(:))>1 %if there are values in the species matrix that are above 1, i.e. measure is Abundance
        measurement = 'Abundance';
        %find abundant species (at least one plot that contains strictly more than 2 individuals)
        is_ab = max(species_matrix,[],1) > ab_thresh; % find the max of each column / boolean vector indicating abundant species
    else %all values are 1 or below
        
        if isequal(species_matrix,round(species_matrix)) %if matrix is round (contains only 0s and 1s) ---> presence
            %measure is presence
            measurement = 'Presence';
            is_ab = sum(species_matrix) > ab_thresh; %find the sum of each column to find species present in at least 3 plots
            
        else %decimal values between 0 and 1 ---> normalised abundance
            measurement = 'Normalised abundance';
            is_ab = sum(species_matrix) > ab_thresh; %find the sum of each column to find species present in at least 3 plots
        end
    end
    
    
    %build species structure:
     
    species.matrix = species_matrix;
    species.plotNames = plotnames;
    species.names = speciesnames;
    species.isAbundant = is_ab; %boolean vector on whether species are abundant
    species.measurement = measurement; %string

end