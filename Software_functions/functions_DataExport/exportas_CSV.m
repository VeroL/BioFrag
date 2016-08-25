% function exportas_CSV(cellarr,csvfile)

% function that accepts any cell array 'cellarr' and writes it to the file 'csvfile' in CSV format
% if 'cellarr' is not a cell array of strings it is converted.
% if csvfile already exists it is overwritten.

function exportas_CSV(cellarr,csvfile)

    %if the cell array does not contain only strings, convert:
    if ~iscellstr(cellarr)
        %convert everything to string
        f = '%9.12f';  %keep 12 digits after decimal points for numbers
        cellarr = cellfun(@(c) (num2str(c,f)), cellarr, 'UniformOutput', false);
    end


    %write CSV file:

    [nrows, ncols] = size(cellarr);

    fid = fopen(csvfile,'w');

    if fid == -1
        disp(['!!! error: exportas_CSV cannot write in "' csvfile '" probably because the file is already open in another programme.']);
    else
        for row=1:nrows
            fprintf(fid,['%s' repmat(',%s',1,ncols-1) '\n'], cellarr{row,:});
        end

        fclose(fid);
    end

end
