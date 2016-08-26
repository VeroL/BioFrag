%Edge Response sequence of operations. function called by a Main (that may vary depending on users)

function Edge_response_ExecSequence(runvar)

    %input runvar (structure) will contain path and file names that depend on users or input data
    %at this stage it contains plugin_info (empty if no plugin, function handle if there is a plugin)
    % and the first logged hello message + tmp log file name

    try 
        
%1) set the software environment (relative folder paths), and output filename conventions:
        
        %--param structure contains naming conventions that do not vary
        param = set_software_parameters;
        
        %--display hello msg and check users have necessary input files to run software:
        qstring = runvar.start_msg;
        yesbutton = 'OK';
        nobutton = 'Cancel';
        button = questdlg(qstring,param.msgbox_Title,yesbutton,nobutton,yesbutton);
        switch button
            case yesbutton %user has files
                %nothing to do
            case nobutton %user wants to terminate run
                errID = 'userterminatedrun';
                assert(false,[param.errorID errID],' '); 
            otherwise %user also wants to terminate run
                errID = 'userterminatedrun';
                assert(false,[param.errorID errID],' ');                
        end
                
       
        %--grab user path (update param)        
        runvar = set_software_environment(runvar, param); %if error here write log in run path)
        
        %--create and write log (in user path)
        runvar = start_log(runvar);


%2) import user data
        runvar = import_prepare_inputData(runvar, param);

        % check if there is a plugin to transform the user input data (plugin name passed on from the main, if any):
        if ~isempty(runvar.import_plugin) % if field import plugin has been specified in the main
            %grab handle
            plugin_handle = runvar.import_plugin;
            %run plugin
            runvar = plugin_handle(runvar, param);    
        end

%3) compute distance map and plot distance values
        msg = '... computation of Binary map (to compute distance to nearest edge)'; dispwrite_log(runvar, param, msg)
        runvar = MakeDistMap(runvar, param);
        msg = '... computation of Distance to nearest Edge for each census point (comparison only)'; dispwrite_log(runvar, param, msg)
        runvar = getPlotDistances(runvar, param);

%4) compute EI map, smooth abundances at plot locations for PCEI values, and landscape metrics
        msg = '... computation of Edge Influence map'; dispwrite_log(runvar, param, msg)
        runvar = MakeEIMapIO(runvar, param);
        msg = '... computation of Edge Influence for each census point + smoothed species abundance'; dispwrite_log(runvar, param, msg)
        runvar = getPlotAbPCEI(runvar, param);
        msg = '... computation of landscape metrics'; dispwrite_log(runvar, param, msg)
        runvar = getLandscape_metrics(runvar, param);

%5) construct training set patterns for PCHD (if file in Matlab folder does not already exists)
        msg = '... computation of Training set of species Edge response'; dispwrite_log(runvar, param, msg)
        runvar = construct_AbPatRef7(runvar, param);

%6) select training set patterns according to experimental design, and rate
        msg = '... computation of assessment of census points distribution wrt PC and EI'; dispwrite_log(runvar, param, msg)
        runvar = make_datasetassessment2(runvar, param);
        %computes plot range in PC/HD by checking if patterns are potential (plots wide apart enough to identify patterns) 
        %creates plotrange and pattern potentiality that species classification needs
        %stores patterns value at plot locations (only for valid patterns)

%7) classify species with closest pattern class, and compute fragmentation impact with PC/HD
        % + generate figures of PCHD graph +FI
        % only analyses/ make figure if species actually has some abundance, not zero everywhere, waste of time and space
        msg = ['... computation of species Edge response and Fragmentation Impact + saving figures.' ...
            ' This operation may take a while but you can start reviewing already saved figures in: ' runvar.folder_Data_Output]; dispwrite_log(runvar, param, msg)
        runvar = classifySpecies_PCHD_dev7(runvar, param);


%8) generate outputs: tif file of EI (+ bin and dist), csv with categories and dataset rating (+ smooth abundance, and plot properties
        msg = '... generating output files tif and csv'; dispwrite_log(runvar, param, msg)
        generate_outputs(runvar, param);

%9) display successful completion message. (messages are written to file as they are generated (using dispwrite_log(param, msg) function)

        msg=[' ****** Edge response completed successfully for ' runvar.PCmap_name ', ' runvar.plot_name ', and ' runvar.species_name]; dispwrite_log(runvar, param, msg)
        msg=['You can view the Edge response outputs in ' runvar.folder_Data_Output '. Bye for now']; dispwrite_log(runvar, param, msg)
    
        
    catch theError  %what to do if an error occurs during run
               
        switch theError.identifier
            case [param.errorID 'inputs'] %if the error is one I generate of type 'inputs' -------------------------------------------------
                
                errmsg = {'There seems to be a problem with the input files selected.', ...
                    'The Edge response software run will terminate.', ...
                    ['Please read the log file ' runvar.log_filename ' for a description of the problem.']};
                
                
                
                uiwait(msgbox(errmsg,param.msgbox_Title,'error','modal')); %wait for user to click ok (uiwait) and force window to stay on top (modal)
                %we don't want to log this one
                
                errmsg = {'If you are unable to solve this issue', ...
                    ['please send the log file: ' runvar.log_filename], ...
                    'along with your input data files (map, coordinates and species matrix)', ...
                    'to Veronique Lefebvre: vero.a.lefebvre@gmail.com', ...
                    'Bye for now'};
                uiwait(msgbox(errmsg,param.msgbox_Title,'error','modal'));
                fid = fopen(runvar.log_filename, 'a');
                nrows= length(errmsg);
                for row=1:nrows
                    fprintf(fid, '%s \r\n', errmsg{row});
                end
                fclose(fid);
                
                
                
                
            case [param.errorID 'userterminatedrun'] %if the error is one I generate of type 'userterminatedrun' ---------------------------
                
                errmsg ={['The ' param.msgbox_Title ' run will terminate. Bye for now']};
                uiwait(msgbox(errmsg,param.msgbox_Title,'error','modal')); %the user wants to terminate the run
                
                
                
                
                
                
            otherwise % any other "unknown" error ------------------------------------------------------------------------------------------
                
                %display and log the error message generated by matlab
                errmsg = ['!!! ' getReport(theError,'extended','hyperlinks','off')];
                uiwait(msgbox(errmsg,param.msgbox_Title,'error','modal')); %wait for user to click ok (uiwait) and force window to stay on top (modal)
                fid = fopen(runvar.log_filename, 'a');
                fprintf(fid, '%s \r\n', errmsg);
                fclose(fid);
                
                %display and log final message
                errmsg = {'An unforeseen error has occured.', ...
                    'The Edge response software run will terminate.', ...
                    ['If you are unable to solve this issue please send the log file: ' runvar.log_filename], ...
                    'along with your input data files (map, coordinates and species matrix)', ...
                    'to Veronique Lefebvre: vero.a.lefebvre@gmail.com', ...
                    'Bye for now'};
                
                uiwait(msgbox(errmsg,param.msgbox_Title,'error','modal'));
                fid = fopen(runvar.log_filename, 'a');
                nrows= length(errmsg);
                for row=1:nrows
                    fprintf(fid, '%s \r\n', errmsg{row});
                end
                fclose(fid);
                
                
        end %end switch on error type
        
        close all;
        
    end %end try catch

end %end function


