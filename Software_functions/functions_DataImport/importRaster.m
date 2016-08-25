% importRaster

% import *geotiff* file (should be a projected map)
% saves map, reference object and tiff info

function [Lmap, R, ginfo, res, ismap_geotif] = importRaster(filename)
    
    %set the last warning to an empty string
    lastwarn('');
    
    %try to get the geo info from the geotif file
    ginfo=geotiffinfo(filename);
    
    %check if Matlab has issued a 'cannot find geo info' warning during the call to geotiffinfo :
    [~, warnid] = lastwarn;
    
    if strcmp(warnid,'map:geotiff:undefinedGTModelTypeGeoKey')
        ismap_geotif = false;
    else
        ismap_geotif = true;
    end
    

    %check if image is indexed 
    %get map and spatial referencing object (MapRasterReference class) - > i.e. projected
    if strcmp(ginfo.ColorType,'indexed')
        [Lmap, ~, R] = geotiffread(filename);
    else
        [Lmap, R] = geotiffread(filename);
    end
    
    %extract resolution
    res = ginfo.PixelScale(1);
    
    %make sure map is double
    Lmap=double(Lmap);
   
       
end