%MakeEIMap
%compute gaussian filter and threshold it at DEI
%compute AVG and SD from the point cover map 
%compute EI map as max deviation

function [HDmap, SDmap, LCmap] = MakeEIMap(PCmap, res, DEI)
       
    %---------------compute thresholded gaussian filter, of DEI size

    %DEI size
    sigma = round(DEI /(2*res)); %gaussian sigma in pixels, half circular radius

    %window size. the default in Matlab2015 imgauusfilt is: 2*ceil(2*SIGMA)+1
    wdsz = 2*ceil(2*sigma)+1; 
    %gaussian filter:    
    h=fspecial('gaussian',wdsz,sigma);
    
    %set the gaussian filter to zero beyond DEI radius, using a circular mask:
    circ_mask = zeros(wdsz);
    circ_mask(ceil(wdsz/2),ceil(wdsz/2))=1; circ_maskaux = bwdist(circ_mask); %compute distance transform from centre
    circ_mask = double(circ_maskaux <= 2*sigma);  %circular mask is 1 within DEI/res pixels from centre, 0 outside
    
    h = h .* circ_mask; % apply mask to guassian filter (threshold)
    h = h / sum(h(:));  %ensure sum of gaussian filter is one

    %-----------------computing spatial average (Local Cover LC)
    
    LCmap = imfilter(PCmap,h,'replicate');   % also tested circular, 0, and symmetric
    
    
    %-----------------computing local standard deviation from E[X^2] - E[X]^2

    h = h*numel(h); %scale h so that sum(h(:)) is not 1 (and fspecial produces filters whose sum of values = 1)
    n = sum(h(:)); n1 = n - 1;

    conv1 = imfilter(PCmap.^2, h/n1 , 'symmetric'); %local average of square image
    conv2 = imfilter(PCmap, h, 'symmetric').^2 / (n*n1); %square local average of image
    SDmap = sqrt(max((conv1 - conv2),0)); % std is sqrt(E[X^2] - E[X]^2)
   
    
    %-----------------computing edge influence (previously called HD for Habitat Disturbance)
    %Point deviation 
    PDmap = LCmap - PCmap;     
    %Edge influence, max of abs(LC-PC) and std, with sign of LC-PC
    HDmap = max(abs(PDmap), SDmap) .* sign(PDmap);
    
end