%% sample the pair of super pixels
close all; clearvars
addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
superpixel_input_dir = fullfile(sourcedir,'super_pixels');

listDir = dir(fullfile(superpixel_input_dir,'*_se1_minNuc3_minStr5_minLum10_objects'));
imagepaths = {listDir.name}';
numImages = length(imagepaths);%


ObjectsMapFileName = imagepaths{2}; 
im_splitStr = regexp(ObjectsMapFileName,'\_','split');
% read the file with nearest neighbor id
objectFeatureFname = fullfile(superpixel_input_dir,[im_splitStr{1} '_superpixel_feature']);
objectFeatures = dlmread(objectFeatureFname, ',',2,0);
features = objectFeatures(:,2);
% read the file with the object feature
nearestNeighborFname = fullfile(superpixel_input_dir,[im_splitStr{1} '_se1_minNuc3_minStr5_minLum10_first50neigh_IDs']);
nearestNeighbor = dlmread(nearestNeighborFname, ',',1,1);
numObjects = size(nearestNeighbor,1);
numNeighbors = 20; % randomly go 20 down the list
sampleSP = zeros(numObjects*numNeighbors,2);
% loop through all the object
for obj = 1:numObjects
    %neigh_indx = randperm(numNeighbors); % get the neighbor index permuted
    %neigh = nearestNeighbor(obj,neigh_indx(1));
    sampleSP(((obj-1)*numNeighbors+ 1):(obj*numNeighbors),1) = objectFeatures(obj,2);
    if sum(nearestNeighbor(obj,:)) == 0
        sampleSP(((obj-1)*numNeighbors + 1):(obj*numNeighbors),1) = objectFeatures(obj,2);
    else
        sampleSP(((obj-1)*numNeighbors + 1):(obj*numNeighbors),2) = features(nearestNeighbor(obj,1:numNeighbors))';
        %sampleSP(obj,2) = objectFeatures(nearestneighbor,2);
    end
end
figure; ndhist(sampleSP(:,1),sampleSP(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 

% normal pixel
opts = setEnvironment_affinity;
opts.features.which_features = {'hue opp'};
which_features = opts.features.which_features;
Nsamples = 10000;%numObjects;

raw_image = imread(fullfile(tiles_dir,[im_splitStr{1} '.tif']));
I = double(raw_image);
f_maps = getFeatures(I,1,which_features,opts);
f_map = f_maps{1};
opts.sigma = 0.5; F1 = sampleF(f_map,Nsamples,opts);
opts.sigma = 1; F2 =  sampleF(f_map,Nsamples,opts);
opts.sigma = 5; F3 =  sampleF(f_map,Nsamples,opts);
F = [F1; F2; F3];
figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 

