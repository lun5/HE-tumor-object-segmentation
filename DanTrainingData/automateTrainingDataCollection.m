% script to extract training data 
% Luong Nguyen 9/14/15

svs_name = 'tp10-867-1';
im_dir = 'Z:\TilesForLabeling';
fileList = dir(fullfile(im_dir,[svs_name '*.tif']));
im_names = {fileList.name}';

Hpixels = [];
Epixels = []; 

for i = 1:length(im_names)
    imname = im_names{i};
    I = imread(fullfile(im_dir,imname));
    [~, H, E] =  extractStainPixels(I,[],0.2,[1 1.15]);
    Hpixels = cat(1,Hpixels, H);
    Epixels = cat(1,Epixels, E);
end