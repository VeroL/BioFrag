

function runvar = start_log(runvar)

    %create log file, and log start message:
    
    %grab current time to create log file name
    tsp = datestr(now, 'yymmdd_HHMMSS');
    runvar.log_filename = fullfile(runvar.folder_Data_Output, ['Edge_response_log_' tsp '.txt']); %overwrite log file name
    
    %write start msg to file:
    fid = fopen(runvar.log_filename, 'w');
    if iscell(runvar.start_msg)
        nrows= length(runvar.start_msg);
        for row=1:nrows
            fprintf(fid, '%s \r\n', runvar.start_msg{row});
        end
    else
        fprintf(fid, '%s \r\n', runvar.start_msg);
    end
    fclose(fid);

end