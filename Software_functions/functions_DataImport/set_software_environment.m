%set_software_environment

%adds the user path to the runvar struct - user path should not normally vary for different runs by the same user.

%runvar contains path and file names that depend on input data
%param contains naming conventions that do not vary between runs

function runvar = set_software_environment(runvar, param)


    %1) grab the user path for the software data

    %grab the location from which the software is run:
    %e.g. if run from a desktop shortcut, or from .exe wihtin the unpacked location
    current_folder = pwd;
    
    %check if the Edge_response_Data_Input folder is located in the current_folder
    
    if  exist(fullfile(current_folder,param.foldername_Data_Input),'dir') %software is run from its folder
        %the software folders are located in the run folder, so proceed
        runvar.software_path = current_folder;
        %inform the user
        msg = ['The Data subfolders are located in : ' runvar.software_path]; 
        uiwait(msgbox(msg,param.msgbox_Title,'help','modal')); %wait for user to click ok (uiwait) and force window to stay on top (modal)
        
    else %software is run from different location --> check with user where they want to store the software data
        
        qstring = {['Do you want the ' param.msgbox_Title ' software to store its Data subfolders in this folder:'];
            [current_folder ' ?']};
        yesbutton = 'Yes, proceed';
        nobutton = 'No, browse';

        button = questdlg(qstring,param.msgbox_Title,yesbutton,nobutton,yesbutton);

        switch button
            case yesbutton %user happy with location
                runvar.software_path = current_folder; 
            case nobutton  %user wants to browse to a different location
                %display folder selection dialog
                software_path = uigetdir(current_folder,[param.msgbox_Title ' : Select a folder to store the ' param.msgbox_Title ' software subfolders:']);
                %prepare to throw error if no folder selected:
                var_to_assert = ischar(software_path);
                erreason = ['NO folder selected to store the ' param.msgbox_Title ' Data subfolders. Log file will be stored in : ' strrep(current_folder, '\', '/')];
                errID = 'inputs';
                check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run if no folder selected
                
                %if no error then set path
                runvar.software_path = software_path;
                
                %inform the user of their path selection, if they made a mistake -> too bad.
                if  ~exist(fullfile(runvar.software_path,param.foldername_Data_Matlab),'dir') 
                    msg = ['The Data subfolders will be stored in your selected folder: ' runvar.software_path]; 
                else
                    msg = ['The Data subfolders are already in the selected folder: ' runvar.software_path]; 
                end
                uiwait(msgbox(msg,param.msgbox_Title,'help','modal')); %wait for user to click ok (uiwait) and force window to stay on top (modal)
                
            otherwise %user clicked cancel - no folder selected
                %throw error
                var_to_assert = 0;
                erreason = ['NO folder selected to store the ' param.msgbox_Title ' Data subfolders. Log file will be stored in : ' strrep(current_folder, '\', '/')];
                errID = 'inputs';
                check_log_error(var_to_assert, erreason, errID, runvar, param); %terminates run as no folder selected
        end
    
    end
    
   
    %2) set paths to working data folders and IO data folders    
    runvar.folder_Data_Input = fullfile(runvar.software_path,param.foldername_Data_Input);
    runvar.folder_Data_Matlab = fullfile(runvar.software_path,param.foldername_Data_Matlab);
    runvar.folder_Data_Output = fullfile(runvar.software_path,param.foldername_Data_Output);
    
    
    %3) create folders when necessary:
    %check if the Data_Input folder exists, and if not create
    if  ~exist(runvar.folder_Data_Input,'dir') %use absolute path to check existence
        mkdir(runvar.folder_Data_Input);       
    end
    if  ~exist(runvar.folder_Data_Matlab,'dir') 
        mkdir(runvar.folder_Data_Matlab);       
    end
    if  ~exist(runvar.folder_Data_Output,'dir') %use absolute path to check existence
        mkdir(runvar.folder_Data_Output);       
    end

end