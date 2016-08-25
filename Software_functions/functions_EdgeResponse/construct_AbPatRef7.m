%making training set: simulated species abundance on PC/EI graph - per category (Forest Core,Forest Edge, etc..)
% using combinations of flipped 1D sigmoids
%these are graphs: the lowest values EI values correspond to the lowest row values so they are at the top
%needs to be displayed inverted (axis xy, not image)


function runvar = construct_AbPatRef7(runvar, param)
   
    runvar.AbPatRef_matfile = fullfile(runvar.folder_Data_Matlab, param.AbPatRef);
    
    if ~exist(runvar.AbPatRef_matfile, 'file') %if file of training set does not exist then compute

        %-------------------------------------------------1) set parameters + make base 1D sigmoids
        
        yEI = -100:100; nbY = numel(yEI);
        xPC = 0:100; nbX = numel(xPC);
        
        propAbThresh = 2/3; %max zone

        amp=[0.25 0.5 1];
        
        fr = 0.25; %fast rate
        mr = 0.15; %medium rate
        sr = 0.1; %slow rate

        %---------PC
        infpointPC    = [55 75 95 65 85]; 
        ratePC        = [fr fr fr mr sr];
        tagPC         = {'Sp','Ds','vDs','Sp','Ds'}; %for sparse dense very_dense
        widthPC = [10 130]; 
        valsPC = xPC;
        [subsigformPC, subsigtagPC, nb_subsigPC] = construct_subsig(infpointPC, ratePC, tagPC, widthPC, valsPC);
        
        %---------Edge
        infpointEdge    = [15 30 45 45 60 75 85]; 
        rateEdge        = [fr mr sr fr mr fr sr];
        tagEdge         = {'T3','T3','T2','T2','T1','T1','T1'};
        widthEdge = [10 130]; 
        valsEdge = yEI;
        [subsigformEdge, subsigtagEdge, nb_subsigEdge] = construct_subsig(infpointEdge, rateEdge, tagEdge, widthEdge, valsEdge);
        
        %---------Core
        infpointCore    = -[5 15 25 10];
        rateCore        = [fr fr fr sr];
        tagCore         = repmat({'0'},[numel(infpointCore) 1]); %0 - ony tag for core
        widthCore =  200; %long width out of spectrum
        valsCore = -abs(yEI);
        [subsigformCore, subsigtagCore, nb_subsigCore] = construct_subsig(infpointCore, rateCore, tagCore, widthCore, valsCore);
        
        %---------Flat
        flatamp = {0.7, [0.5 0.6], [0 0.05 0.1], 0}; %generalist noPref, unknown, nonabundant, absent
        
         
        %-------------------------------------------------2) make PC and EI 1D sigmoid bank + corresponding category name
        
        %---------PC
        sigformPC = [subsigformPC;      %Forest
            fliplr(subsigformPC)];      %Matrix (flip on x axis to get Matrix) 
        signamePC = [repmat({'Forest'},[nb_subsigPC,1]);    %adding name to each sigform
            repmat({'Matrix'},[nb_subsigPC,1])]; 
        sigtagPC = [subsigtagPC;                            %adding tag to each sigform
            subsigtagPC];
                
        
        %---------EI
        sigformEI = [subsigformCore;    %Core
            subsigformEdge;                %Edge lowside
            fliplr(subsigformEdge);     %Edge highside
            ones(1,nbY)];               %noPref
        signameEI = [repmat({'Core'},[nb_subsigCore,1]);    %adding name to each sigform 
            repmat({'Edge'},[nb_subsigEdge,1]);
            repmat({'Edge'},[nb_subsigEdge,1]);
            {'noPref'}]; 
        sigtagEI = [subsigtagCore;                            %adding tag to each sigform
            cellfun(@(x) {[x ' lowside']}, subsigtagEdge);
            cellfun(@(x) {[x ' highside']}, subsigtagEdge);
            {'inf'}];

        
        %-------------------------------------------------3) make EI GEN 1D sigmoid bank + corresponding category name
        
        %GEN - EI only
        %no PC filtering
       
        sigformgenEI = [
            0.5*subsigformEdge + fliplr(subsigformEdge)         %forest (highside) edge overflow, half the amplitude (only one combination)
            subsigformEdge + 0.5*fliplr(subsigformEdge)         %other side, matrix (lowside) edge overflow
            subsigformCore(1,:)                                 %only the narrowest core (core generalist)
            subsigformEdge + fliplr(subsigformEdge)             %both sides (edge generalist)
            ]; 
        signamegenEI = [                                        %adding name to each sigform
            repmat({'Forest Edge'},[nb_subsigEdge,1])
            repmat({'Matrix Edge'},[nb_subsigEdge,1])
            {'Generalist Core'} 
            repmat({'Generalist Edge'},[nb_subsigEdge,1])
            ];
        sigtaggenEI = [
            cellfun(@(x) {[x ' overflow']}, subsigtagEdge);    %adding tag to each sigform
            cellfun(@(x) {[x ' overflow']}, subsigtagEdge);
            subsigtagCore(1);
            subsigtagEdge];
        
        
        %FLATS - (generalist noPrefs, unknowns ... )
        sigformflat = [
            flatamp{1}'*ones(1,nbY)                             %matrix product, flatamp rows and nbY cols
            flatamp{2}'*ones(1,nbY)
            flatamp{3}'*ones(1,nbY)
            flatamp{4}'*ones(1,nbY)
            ];

        signameflat = [
            repmat({'Generalist noPref'},[numel(flatamp{1}),1])
            repmat({'Unknown'},[numel(flatamp{2}),1])
            repmat({'Nonabundant'},[numel(flatamp{3}),1])
            repmat({'Absent'},[numel(flatamp{4}),1])
            ];

        
        %-------------------------------------------------4) combine 1D sigmoid to make image pattern and store in structure
        
        %number of sigmoids
        nb_sigPC = size(sigformPC,1);
        nb_sigEI = size(sigformEI,1);
        nb_siggenEI = size(sigformgenEI,1);
        nb_flat = size(sigformflat,1);
               
        nbpatterns= nb_sigPC*nb_sigEI + nb_siggenEI + nb_flat; 
        
        %memory allocation:

        AbPatRef.Habitat = [];            %Habitat name
        AbPatRef.EIPref = [];             %Edge Influence Pref 
        AbPatRef.Category = [];           %category name composed of Habitat and Edge effects pref 
        AbPatRef.ScaleCategory = [];      %category name composed of Habitat and Edge effects pref tags and amplitude. used in naive bayes
        AbPatRef.isScalemax = [];         %1 if pattern has the max scale of category (amp = 1)
        AbPatRef.matrix = [];             %values of modelled abundance for all PC/EI pairs (matrix)
        AbPatRef.matrix_maxbin = [];      % 1 if abundance value is over proportion threshold propAbThresh ( > 2/3 of amp), 0 otherwise
        %copy over 
        AbPatRef(2:nbpatterns,:)=AbPatRef(1);
        indpat = 0;

        %2D graphs indices
        [pcvals, eivals] = meshgrid(xPC, yEI);
        
        %mask of possible combinations of PC and EI
        graphmask = eivals<110-pcvals & eivals>-pcvals-10;
        
        %------------------------------------make patterns by combining EI and PC - NON GEN patterns
        
        for indPC = 1:nb_sigPC
            for indEI = 1:nb_sigEI

                abpat = repmat(sigformPC(indPC,:),nbY,1) .* repmat(sigformEI(indEI,:)',1,nbX);
                
                %exclude patterns where max abundance is outside possible combination zone of PC and EI:
                if max( max( abpat.* graphmask )) >=0.9  %only keep if max masked pattern above 1 (all patterns are normalised)

                    %loop over amplitude values
                    if strcmp(signameEI{indEI},'noPref')
                        ampvals = amp(2:end); %only used highest amp values for noPref species
                    else
                        ampvals = amp;
                    end
                    for indamp = 1:numel(ampvals)

                        %mulitply normalised pattern by amplitude
                        abpatamp =  abpat*ampvals(indamp);

                        %insert pattern in AbpatRef

                        indpat = indpat + 1;

                        AbPatRef(indpat).Habitat = [signamePC{indPC} ' ' sigtagPC{indPC}];                %Habitat name and tag
                        AbPatRef(indpat).EIPref = [signameEI{indEI} ' ' sigtagEI{indEI}];                 %EI pref name and tag
                        AbPatRef(indpat).Category = [signamePC{indPC} ' ' signameEI{indEI}];              %Main category: habitat and EIpref names 
                        
                        %Naive Bayes (classification) category: habitat and EIpref names and tags + amplitude 
                        AbPatRef(indpat).ScaleCategory = [AbPatRef(indpat).Habitat ' ' AbPatRef(indpat).EIPref ' ' num2str(ampvals(indamp))];             

                        AbPatRef(indpat).isScalemax = ampvals(indamp) == 1;                                     %1 if pattern has the max scale of category
                        AbPatRef(indpat).matrix = abpatamp;                                                     %values of modelled abundance for all PC/EI pairs (matrix)
                        AbPatRef(indpat).matrix_maxbin = bwlabel(abpatamp >= propAbThresh * max(abpatamp(:)));  %label binary matrix ( 1 where abundance is over 2/3 of max abundance)
                                                                                                                %only 1 max zone except for generalists core and edge (this is to make sure there are plots in both zones)
                    end %end loop on amplitude
                end %end condition on mask
                
            end %end double for loop on PC and EI sigmoids
        end
  
        
        %------------------------------------make patterns by combining EI and PC - GEN patterns
        
        %loop on EI sigmoids, no PC filtering
        for indEI = 1:nb_siggenEI

            abpat = repmat(sigformgenEI(indEI,:)',1,nbX);

            %loop over amplitude values
            ampvals = amp(2:end); %only used highest amp values for generalist species

            for indamp = 1:numel(ampvals)

                %mulitply normalised pattern by amplitude
                abpatamp =  abpat*ampvals(indamp);

                %insert pattern in AbpatRef

                indpat = indpat + 1;
                
                patcats = textscan(signamegenEI{indEI}, '%s %s');

                AbPatRef(indpat).Habitat = [patcats{1}{:}];                %Habitat name and tag
                AbPatRef(indpat).EIPref = [patcats{2}{:} ' ' sigtaggenEI{indEI}];                 %EI pref name and tag
                AbPatRef(indpat).Category = signamegenEI{indEI};              %Main category: habitat and EIpref names 

                %Naive Bayes (classification) category: habitat and EIpref names and tags + amplitude 
                AbPatRef(indpat).ScaleCategory = [AbPatRef(indpat).Habitat ' ' AbPatRef(indpat).EIPref ' ' num2str(ampvals(indamp))];             

                AbPatRef(indpat).isScalemax = ampvals(indamp) == 1;                                     %1 if pattern has the max scale of category
                AbPatRef(indpat).matrix = abpatamp;                                                     %values of modelled abundance for all PC/EI pairs (matrix)
                AbPatRef(indpat).matrix_maxbin = bwlabel(abpatamp >= propAbThresh * max(abpatamp(:)));  %label binary matrix ( 1 where abundance is over 2/3 of max abundance)
                                                                                                        %only 1 max zone except for generalists core and edge (this is to make sure there are plots in both zones)
            end %end loop on amplitude                

        end %end for loop on EI gen sigmoids
        
        
        %-------------------------------------------adding flat patterns at the end
        for indsig = 1:size(sigformflat,1)
            abpat = repmat(sigformflat(indsig,:)',1,nbX);
            indpat = indpat + 1;

            patname = signameflat{indsig};
            patcats = textscan(patname, '%s %s');  
            AbPatRef(indpat).Habitat = patcats{1}{:};                                 
            if isempty(patcats{2})
                AbPatRef(indpat).EIPref = patcats{1}{:};                          
            else
                AbPatRef(indpat).EIPref = patcats{2}{:};
            end
            AbPatRef(indpat).Category = patname;  
            AbPatRef(indpat).ScaleCategory = patname;   
            AbPatRef(indpat).isScalemax = max(abpat(:)) == max(flatamp{1});                                                           
            AbPatRef(indpat).matrix = abpat;                                                   
            AbPatRef(indpat).matrix_maxbin = bwlabel(abpat >= propAbThresh * max(abpat(:)));  
        end
        
        
        %---------------------------------- Trim AbPatRef structure
        AbPatRef(indpat+1:end) = []; %remove empty structure elements
        nbpat = length(AbPatRef); %get final number of patterns



        %-------------------------------------------------5)--store general properties across all patterns:

        %number of unique pattern shapes (ignoring intensity scaling) - used to assess plot spread :
        AbPatRef_properties.nbuniqzones = sum([AbPatRef.isScalemax]);

        %and their index in AbPatRef :
        AbPatRef_properties.ind_uniquezone_pattern = find([AbPatRef.isScalemax])';

        %and index of unique shapes as a boolean vector (1 if unique / 0 if duplicate shape
        AbPatRef_properties.bool_uniquezone_pattern = [AbPatRef.isScalemax]';

        %store list of main categories:
        AbPatRef_properties.main_category_list = unique({AbPatRef.Category}','stable'); 
        %     AbPatRef_properties.main_category_list = {'Forest Core','Forest Edge','Forest noPref','Matrix Core','Matrix Edge','Matrix noPref', ...
        %      'Generalist Core','Generalist Edge','Generalist noPref','Unknown','Nonabundant'};

        %store list of sub categories:
        subcats = arrayfun(@(x) [x.Habitat ' ' x.EIPref], AbPatRef, 'UniformOutput', false)';
        AbPatRef_properties.sub_category_list = unique(subcats','stable'); 

        %store sum of abundance for each pattern
        AbPatRef_properties.sum_abpat = zeros(nbpat,1);
        for indpat = 1:nbpat
            AbPatRef_properties.sum_abpat(indpat) = sum(AbPatRef(indpat).matrix(:));
        end
        
        
        
        
        %------------------------------------------------------save structure array  and its properties as .mat file
        save(runvar.AbPatRef_matfile,'AbPatRef','AbPatRef_properties');
        msg=['Training set of species Edge response was computed ' runvar.AbPatRef_matfile]; dispwrite_log(runvar, param, msg)
        
    else %if exists just load:
        msg=['previously computed: Training set of species Edge response ' runvar.AbPatRef_matfile]; dispwrite_log(runvar, param, msg)
    end %end of need_compute condition   
    
    %if need_fig, fig_AbPatRef(AbPatRef); end %make fig if needed
    

end %end of construct_AbPatRef function




%---------------------------- function to create the sigmoid forms

function [subsigform, subsigtag, nb_subsig] = construct_subsig(infpoint, rate, tag, width, vals)

    %infpoint, rate and tag should be same length and correspond
    ni = numel(infpoint);
    nw = numel(width);
    
    %store 1D sigmoids
    nb_subsig = ni * nw;
    subsigform = zeros(nb_subsig, numel(vals));
    subsigtag = cell(nb_subsig,1);
    for wi = 1:nw
        for ii=1:ni
            indsig = (wi-1)*ni + ii;
            subsigform(indsig,:) = 1./(1+exp( -rate(ii) * (vals - infpoint(ii) ) )) ...
                    - 1./(1+exp( -rate(ii) * (vals - (infpoint(ii) + width(wi)) ) )); 
            subsigform(indsig,:) = subsigform(indsig,:)/max(subsigform(indsig,:));
            
            subsigtag{indsig}=tag{ii};
        end
    end
    
end



% %---------------------------- function to display the training set by categories
% 
% function fig_AbPatRef(AbPatRef,subsigformCore,subsigformEdge,subsigtagEdge,subsigformPC,subsigtagPC)
% 
%     xPC = 0:100; 
%     yEI = -100:100;
%     
%     
%     %figures of 1D sigmoids:
%     figure('position',[10 200 1200 400]); plot(yEI,subsigformCore);set(gca,'XTick',-100:10:100);grid on;
%     
%     nbe = size(subsigformEdge,1);
%     figure('position',[10 200 1200 400]); plot(yEI,subsigformEdge(1:nbe/2,:) );set(gca,'XTick',-100:10:100);grid on; legend(subsigtagEdge(1:nbe/2),'Location','bestoutside');
%     figure('position',[10 200 1200 400]); plot(yEI,subsigformEdge(nbe/2+1:nbe,:) );set(gca,'XTick',-100:10:100);grid on; legend(subsigtagEdge(nbe/2+1:nbe),'Location','bestoutside');
% 
%     nbp = size(subsigformPC,1);
%     figure('position',[10 200 1200 400]); plot(xPC,subsigformPC(1:nbp/2,:) );set(gca,'XTick',0:10:100);grid on; legend(subsigtagPC(1:nbp/2),'Location','bestoutside');
%     figure('position',[10 200 1200 400]); plot(xPC,subsigformPC(nbp/2+1:nbp,:) );set(gca,'XTick',0:10:100);grid on; legend(subsigtagPC(nbp/2+1:nbp),'Location','bestoutside');
% 
%     
%     %figures of 2D sigmoids
%     [pcvals, eivals] = meshgrid(xPC, yEI);
% 
%     [catlist, ~, catlab] = unique({AbPatRef.Category}','stable');
%     
%     for indcat = 1 : length(catlist)
%         indpat_incat_scmax = find( catlab == indcat & [AbPatRef.isScalemax]'==1);
%         nbpat_todisplay = numel(indpat_incat_scmax);
%         numpat=0;
%         while numpat < nbpat_todisplay
%             figure;
%             indpat1 = numpat +1;
%             indpatend = min(numpat+42,nbpat_todisplay);
%             indsub=0;
%             for indpat = indpat1 : indpatend; 
%                 numpat = numpat + 1;
%                 indsub = indsub + 1;
%                 abpatmasked = AbPatRef(indpat_incat_scmax(indpat)).matrix; %_maxbin;
%                 abpatmasked(eivals>110-pcvals | eivals<-pcvals-10)=0;
%                 subplot(6,7,indsub);
%                 imagesc(0:100, -100:100, abpatmasked); colormap jet; axis xy; axis image
%                 hold on; plot([0 100],[0 0],'k'); plot([0 100],[50 0],'k--'); plot([0 100],[0 -50],'k--');
%                 set(gca,'XTick',[]); set(gca,'YTick',[]);
%                 title(AbPatRef(indpat_incat_scmax(indpat)).ScaleCategory,'fontsize',8);
%             end
%         end
%     end
% 
% end
