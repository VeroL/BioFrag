
function dispwrite_log(runvar, param, msg)
 
    msgbox(msg,param.msgbox_Title,'help','replace'); %display box

    %append msg to file:
    fid = fopen(runvar.log_filename, 'a');
    if iscell(msg)
        nrows= length(msg);
        for row=1:nrows
            fprintf(fid, '%s \r\n', msg{row});
        end
    else
        fprintf(fid, '%s \r\n', msg);
    end
    fclose(fid);
end