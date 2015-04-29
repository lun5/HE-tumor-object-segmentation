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
    fname_features = fullfile(superpixel_input_dir,[im_splitStr{1},'_superpixel_feature']);
    features_matrix = dlmread(fname_features,',',1,0);
    features = features_matrix(:,2)';
    fname_neighborIDs = fullfile(superpixel_input_dir,[im_splitStr{1},'_se1_minNuc3_minStr5_minLum10_first50neigh_IDs']);
    neighborIDs = dlmread(fname_neighborIDs,',');
    
    numObjects = size(neighborIDs,1);
    neighbor_features = zeros(size(neighborIDs));
    for indx = 1:numObjects
        neighbors = neighborIDs(indx,:) + 1;
        if sum(neighbors == 0) > 0
            neighbor_features(indx,:) = features(1);
        else
            neighbor_features(indx,:) = features(neighbors);
        end
    end
    
    fname = fullfile(superpixel_input_dir,[im_splitStr{1},'_se1_minNuc3_minStr5_minLum10_first50neigh_huefeatures']);
    dlmwrite(fname,neighbor_features);
    sprintf('Finish with file %s\n',im_splitStr{1}) 
    %toc
end

