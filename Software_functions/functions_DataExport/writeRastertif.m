%write geotiff

function writeRastertif(pathfilename,map,R,ginfo)

    p=ginfo.GeoTIFFTags.GeoKeyDirectoryTag;
    geotiffwrite(pathfilename,map,R,'GeoKeyDirectoryTag',p);
end