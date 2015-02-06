%% script to convert tiles to jpg
% run the code in parallel
%pool = gcp;
addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling');
jpeg_dir = fullfile(sourcedir,'TilesForLabeling_jpg');

if ~exist(jpeg_dir,'dir')
    mkdir(jpeg_dir);
    fileattrib(jpeg_dir,'+w');
end

fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths); 
for i = 1: numImages
    imname = imagepaths{i};
    splitStr = regexp(imname,'\.','split');
    raw_image = imread(fullfile(tiles_dir,imname));
    imwrite(raw_image,fullfile(jpeg_dir,[splitStr{1},'.jpg']),'jp2','Mode','lossless');
end

