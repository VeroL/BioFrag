%MakeDistMap 
%compute dist map from the PC map

function runvar = MakeDistMap(runvar, param)

    %create filenames for binary map and dist map
    runvar.Binmap_matfile = fullfile(runvar.folder_Data_Matlab, [runvar.PCmap_name param.PCBin]);
    runvar.Distmap_matfile = fullfile(runvar.folder_Data_Matlab, [runvar.PCmap_name param.Distance2nE]);
   
    
    %1)--------create or load binary map
    if ~exist(runvar.Binmap_matfile,'file')
        
        %load matlab point cover map
        load(runvar.PCmap_matfile);
        
        if isempty(runvar.PCT) %compute threshold
    
            mi=min(PC.map(:)); ma=max(PC.map(:));
            tetI=PC.map-min(PC.map(:));
            tetI=tetI/max(tetI(:));
            V2 = graythresh(tetI);
            V2=V2*(ma-mi)+mi;
            PCBin.threshold = round(V2); %that's the computed threshold
        else %store threhold in var
            PCBin.threshold = runvar.PCT;
        end
        
        %threshold map
        PCBin.map=PC.map>PCBin.threshold;
        
        PCBin.R = PC.R;
        PCBin.ginfo = PC.ginfo;
        PCBin.res = PC.res; 
        %save PCBin structure
        save(runvar.Binmap_matfile,'PCBin');
        
        msg=['computed binary map with threshold ' num2str(PCBin.threshold) ' for map ' runvar.PCmap_name]; dispwrite_log(runvar, param, msg)
         
    else
        load(runvar.Binmap_matfile);
        msg=['previously computed: Binmap ' runvar.Binmap_matfile]; dispwrite_log(runvar, param, msg)
    end
    
    
    msg = '... computation of Distance to nearest Edge map (comparison only)'; dispwrite_log(runvar, param, msg)
        
                
    %2)--------compute distance map from bin map or load
    if ~exist(runvar.Distmap_matfile,'file')
    
        %perforation diameter (to compute distance transform)
        perf_diameter = 60; %in m

        res = PCBin.res;

        %remove perforation smaller than ~60m wide from the clipped binary map
        perf_w = floor(perf_diameter/res); %perforation max width in pixels (better to floor it)
        SEp = strel('square', perf_w);
        fd = imdilate(PCBin.map, SEp); %dilate forest
        fc = imreconstruct(imcomplement(fd), imcomplement(PCBin.map),4); %remove outside dilations, background has a connectivity of 4
        Binmap_np = imcomplement(fc); 

        %distance from all forest pixels to nearest non forest (in meters)
        f_distedge = bwdist(~Binmap_np) * res;

        %distance from all non-forest pixels to nearest forest (in meters (negative))
        nf_distedge = bwdist(Binmap_np) * -res;

        %in and out distance to edge
        dist_map = nf_distedge+f_distedge;


        %2------------save distance map

        %distance map
        D2nE.map = dist_map;
        D2nE.perf_diameter = perf_diameter; %#ok<STRNU>

        save(runvar.Distmap_matfile,'D2nE');

        msg = ['computed distance map for map ' runvar.PCmap_name]; dispwrite_log(runvar, param, msg)
    else
        msg = ['previously computed: Dist map ' runvar.Distmap_matfile]; dispwrite_log(runvar, param, msg)
    end
    
end


