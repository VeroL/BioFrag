%MakeEIMap 
%compute Local Cover from the point cover map 
%compute optimal DEI from plot locations and distance to edge 
%compute the EI map

function runvar = MakeEIMap(runvar, param)

%---------------get the DEI value, set runvar.DEI

    DEIchoices = 250:250:2000;

    %) check whether it has already be given in the main
    if isempty(runvar.DEI)

        %) compute optimum value using distance of plots
        load(runvar.plot_matfile);
        %distance to edge of all plots   
        plot_dist = [plotLoc.dist];
        % ------- compute optimal DEI using distance to nearest edge:
        %this optimum is computed from the distance.. which may not be accurate...
        nbf = 5; % nb of furthest away plots 
        plot_dist_sort = sort(plot_dist); 
        Distopt = round(mean(plot_dist_sort(end-nbf+1:end))); %avg distance of the 5 plots furthest away in forest (dist positive)
        [~, indDEIopt] = min(abs(DEIchoices - Distopt));

        DEIopt = DEIchoices(indDEIopt);      

        %) ask user if automatically selected value is ok (to get core forest species)
        qstring = {sprintf('The optimum depth of edge influence (DEI) so that at least some of your census points are in the forest core area while retaining the deepest possible edge zone is %dm.',Distopt); ...
        sprintf('We suggest using DEI = %dm to compute the Edge Influence map.', DEIopt);
        sprintf('Do you wish to use DEI = %dm ?', DEIopt)};
        yesbutton = sprintf('OK, use DEI = %dm',DEIopt);
        nobutton = 'No, select another DEI value';
        button = questdlg(qstring,param.msgbox_Title,yesbutton,nobutton,yesbutton);
        switch button
            case yesbutton %user is ok with DEIopt
                runvar.DEI = DEIopt;
                msg=sprintf('A Depth of Edge Influence (DEI) of %dm has been selected to match the experimental design.',runvar.DEI); dispwrite_log(runvar, param, msg);
            case nobutton %-------------------------------------user wants to select another value
                %menu
                qstring = {'Select a Depth of Edge Influence value for the computation of the Edge Influence map:                                                   ';
                    sprintf('DEI > %dm => your census points are unlikely to correspond to forest core area, species are unlikely to be classified as "forest core".', DEIopt);
                    sprintf('DEI < %dm => the width of the edge zone is reduced.                                                                                    ', DEIopt)};
                menuchoices = cellstr(num2str(DEIchoices'));
                DEIselected = menu(qstring,menuchoices);
                if DEIselected > 0
                    runvar.DEI = DEIchoices(DEIselected);
                    msg=sprintf('You have selected a Depth of Edge Influence (DEI) of %dm.',runvar.DEI); dispwrite_log(runvar, param, msg);
                else %user did not answer the question
                    errID = 'userterminatedrun';
                    assert(false,[param.errorID errID],' ');
                end
            otherwise %user did not answer the question
                errID = 'userterminatedrun';
                assert(false,[param.errorID errID],' ');                
        end
    else
        msg=sprintf('The Depth of Edge Influence (DEI) used to compute the map of Edge Influence is %dm.',runvar.DEI); dispwrite_log(runvar, param, msg);
    end
         
    
    %DEI currently expected to be a scalar:
    DEIvec = runvar.DEI;
    


    %create filename for Local Cover, SD and EI maps
    runvar.LCmap_matfile = fullfile(runvar.folder_Data_Matlab, [runvar.PCmap_name '_DEI_' num2str(DEIvec) param.LocalCover]);
    
    %-------------------------------compute LC, SD, EI maps from PC map, or load
    if ~exist(runvar.LCmap_matfile,'file')

        %load matlab point cover map
        load(runvar.PCmap_matfile);
        
        %Point Cover (PC):    
        PCmap = PC.map;

        %resolution
        res = PC.res;     

        % -------------------------------------- compute local covers 

        nbDEI = numel(DEIvec);

        %preallocate memory
        LC.map = zeros(size(PCmap));
        LC.DEI = 0;
        LC(2:nbDEI)=LC(1); %copy over 

        SD.map = zeros(size(PCmap));
        SD.DEI = 0;
        SD(2:nbDEI)=SD(1); %copy over 

        HD.map = zeros(size(PCmap));
        HD.DEI = 0;
        HD(2:nbDEI)=HD(1); %copy over 

        for k = 1:nbDEI
            
            %---------------compute thresholded gaussian filter, of DEI size

            %DEI size
            sigma = round(DEIvec(k) /(2*res)); %gaussian sigma in pixels, half circular radius

            %window size. the default in Matlab2015 imgaussfilt is: 2*ceil(2*SIGMA)+1
            wdsz = 2*ceil(2*sigma)+1; 
            %gaussian filter:    
            h=fspecial('gaussian',wdsz,sigma);
            
            %set the gaussian filter to zero beyond DEI radius, using a circular mask:
            circ_mask = zeros(wdsz);
            circ_mask(ceil(wdsz/2),ceil(wdsz/2))=1; circ_maskaux = bwdist(circ_mask); %compute distance transform from centre
            circ_mask = double(circ_maskaux <= 2*sigma); %circular mask is 1 within DEI/res pixels from centre, 0 outside
            
            h = h .* circ_mask; % apply mask to guassian filter (threshold)
            h = h / sum(h(:)); %ensure sum of gaussian filter is one
            
            %-----------------computing spatial average (Local Cover LC)
            
            LC(k).map = imfilter(PCmap,h,'replicate');   % also tested circular, 0, and symmetric
            %store DEI used to produce LC map
            LC(k).DEI = DEIvec(k);
        
            %-----------------computing local standard deviations from E[X^2] - E[X]^2
        
            h = h*numel(h); %scale h so that sum(h(:)) is not 1 (and fspecial produces filters whose sum of values = 1
            n = sum(h(:)); n1 = n - 1;

            conv1 = imfilter(PCmap.^2, h/n1 , 'symmetric'); %local average of square image
            conv2 = imfilter(PCmap, h, 'symmetric').^2 / (n*n1); %square local average of image
            SD(k).map = sqrt(max((conv1 - conv2),0)); % std is sqrt(E[X^2] - E[X]^2)
            %store DEI used to produce SD map
            SD(k).DEI = DEIvec(k);
            
            %-----------------computing edge influence (previously called HD for Habitat Disturbance)
            
            %Point deviation 
            PD = LC(k).map - PCmap;     
            %Edge influence, max of abs(LC-PC) and std, with sign of LC-PC
            HD(k).map = max(abs(PD), SD(k).map) .* sign(PD);
            %store DEI used to produce HD map
            HD(k).DEI = DEIvec(k);

            
        end
        
    
        %------------save LC, SD and EI maps for all DEIs

        save(runvar.LCmap_matfile,'LC','SD','HD');

        %display processing step:
        msg=['computed LC, SD, EI maps with Depth of Edge Influence = ' num2str(DEIvec) ' for map ' runvar.PCmap_name];
        dispwrite_log(runvar, param, msg)
    else
        msg=['previously computed: LC, SD, EI maps with Depth of Edge Influence = ' num2str(DEIvec) ' for map ' runvar.PCmap_name];
        dispwrite_log(runvar, param, msg)
    end
    
end


