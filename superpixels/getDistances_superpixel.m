% calculate the distance from the object to the neighbors
% we can read the X and Y coordinates. Actually, we can read the whole
% matrix, substract that from the object coordinates

%% Calculate the feature of super pixels using the map
% INPUT: f_maps (whatever you give it)
%       _se1_minNuc3_minStr5_minLum10_objects: object ID
% OUTPUT text file with 4 columns (object id, feature values, num pixels) with the name
% of the feature being calculated (for today, let's just hue it)

close all; clearvars
addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
superpixel_input_dir = fullfile(sourcedir,'super_pixels');

listDir = dir(fullfile(superpixel_input_dir,'*_se1_minNuc3_minStr5_minLum10_objects'));
imagepaths = {listDir.name}';
numImages = length(imagepaths);% 420
parfor j = 1:numImages
    %tic;
    ObjectsMapFileName = imagepaths{j}; 
    im_splitStr = regexp(ObjectsMapFileName,'\_','split');
    fname_xcoords = fullfile(superpixel_input_dir,[im_splitStr{1} '_se1_minNuc3_minStr5_minLum10_first50neigh_Xcoors']);
    fname_ycoords = fullfile(superpixel_input_dir,[im_splitStr{1} '_se1_minNuc3_minStr5_minLum10_first50neigh_Ycoors']);
    
    xcoords = dlmread(fname_xcoords,',');
    ycoords = dlmread(fname_ycoords,',');
    
    numNeighs = size(xcoords,2) - 1;
    xdistances = xcoords(:,2:end) - repmat(xcoords(:,1),[1 numNeighs]);
    ydistances = ycoords(:,2:end) - repmat(ycoords(:,1),[1 numNeighs]);
    dists = (xdistances.^2 + ydistances.^2).^.5;
    dists = cat(2,zeros(size(xcoords,1),1),dists);
    fname = fullfile(superpixel_input_dir,[im_splitStr{1},'_se1_minNuc3_minStr5_minLum10_first50neigh_distances']);
    dlmwrite(fname,dists);
    sprintf('Finish with file %s\n',im_splitStr{1}) 
    %toc
end

