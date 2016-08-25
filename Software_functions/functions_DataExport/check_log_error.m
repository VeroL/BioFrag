%prepare error message, and if error has occured display message and log it
%check for error and throw it if needed with ID [param.errorID errID] and message: errmsg

function check_log_error(var_to_assert, erreason, errID, runvar, param)

    %prepare error message:
    errmsg = [param.runaborted erreason];
    
    %if there is an error display a message then log it
    if ~var_to_assert %if there is an error
        uiwait(msgbox(errmsg,param.msgbox_Title,'error','modal')); %wait for user to click ok (uiwait) and force window to stay on top (modal)
        
        %append errmsg to file:
        fid = fopen(runvar.log_filename, 'a');
        if iscell(errmsg)
            nrows= length(errmsg);
            for row=1:nrows
                fprintf(fid, '%s \r\n', errmsg{row});
            end
        else
            fprintf(fid, '%s \r\n', errmsg);
        end
        fclose(fid);
    end
        
    %check for error and throw it if needed with ID [param.errorID errID] and message: errmsg
    assert(logical(var_to_assert),[param.errorID errID],errmsg); 


end