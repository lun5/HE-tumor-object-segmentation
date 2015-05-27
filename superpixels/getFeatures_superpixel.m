%% Calculate the feature of super pixels using the map
% INPUT: f_maps (whatever you give it)
%       _se1_minNuc3_minStr5_minLum10_objects: object ID
% OUTPUT text file with 4 columns (object id, feature values, num pixels) with the name
% of the feature being calculated (for today, let's just hue it)

close all; clearvars
addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
superpixel_input_dir = fullfile(sourcedir,'super_pixels_1');

listDir = dir(fullfile(superpixel_input_dir,'*_se1_minNuc3_minStr5_minLum10_objects'));
imagepaths = {listDir.name}';
numImages = length(imagepaths);% 420
opts = setEnvironment_affinity;
opts.features.which_features = {'hue opp'};
which_features = opts.features.which_features;
tic;
parfor j = 1:numImages
    %tic;
    ObjectsMapFileName = imagepaths{j}; 
    im_splitStr = regexp(ObjectsMapFileName,'\_','split');
    raw_image = imread(fullfile(tiles_dir,[im_splitStr{1} '.tif']));
    I = double(raw_image);
    f_maps = getFeatures(I,1,which_features,opts);
    f_map = f_maps{1};
    objectMap = dlmread(ObjectsMapFileName, ',', 1, 1);
    numObjects = max(objectMap(:));
    superpixel_features = zeros(numObjects + 2,3);
    superpixel_features(1,1) = numObjects;
    superpixel_features(2:end,1) = 0:numObjects;
        
    for indx = 0:numObjects
        features = f_map(objectMap == indx);
        superpixel_features(indx+2,2) = circ_mean(features);
        superpixel_features(indx+2,3) = numel(features); %radius of object
    end    
%     h=figure;
%         subplot(2,1,1);
%         histogram(f_map(:),'Normalization','probability','FaceColor',[0.8 0.8 0.8],'BinWidth',0.2);
%         title('Raw image');set(gca,'FontSize',18);xlim([-pi pi]);   
%         subplot(2,1,2);
%         histogram(superpixel_features(2:end,2),'Normalization','probability','FaceColor',[0.8 0.8 0.8],'BinWidth',0.2);
%         title('Super pixel');set(gca,'FontSize',18);xlim([-pi pi]);                
%         set(gcf,'color','white') % White background for the figure.
%         set(gca,'LooseInset',get(gca,'TightInset'))        

    fname = fullfile(superpixel_input_dir,[im_splitStr{1},'_superpixel_feature']);
    dlmwrite(fname,superpixel_features);
    sprintf('Finish with file %s\n',im_splitStr{1}) 
    %toc
end

