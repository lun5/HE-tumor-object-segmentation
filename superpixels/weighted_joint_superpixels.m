% calculate the joint distribution for pair of pixels
% okay, don't try to make it perfect yet. First try to make it run
close all; clearvars
addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
superpixel_input_dir = fullfile(sourcedir,'super_pixels');

listDir = dir(fullfile(superpixel_input_dir,'*_se1_minNuc3_minStr5_minLum10_objects'));
imagepaths = {listDir.name}';
numImages = length(imagepaths);% 
sigma = 10; % weigh the neighbors by distance
numBins = 50;
binEdges = -pi:0.1:pi; %linspace(-pi,pi,numBins);
binCenters = binEdges + 0.05;
% read the superpixel features break the superpixels into these bins
numNeighbors = 15;
if numNeighbors > 50
    numNeighbors = 50;
end

for j = 1:numImages
    %tic;
    ObjectsMapFileName = imagepaths{j}; 
    im_splitStr = regexp(ObjectsMapFileName,'\_','split');
    
    fname_features = fullfile(superpixel_input_dir,[im_splitStr{1},'_superpixel_feature']);
    features_matrix = dlmread(fname_features,',',2,0);

    fname_neighborFeatures = fullfile(superpixel_input_dir,[im_splitStr{1} ...
        '_se1_minNuc3_minStr5_minLum10_first50neigh_huefeatures']);
    neighborFeatures = dlmread(fname_neighborFeatures,',',1,0);
    
    fname_neighborDistances = fullfile(superpixel_input_dir,[im_splitStr{1} ...
        '_se1_minNuc3_minStr5_minLum10_first50neigh_distances']);
    neighborDistances = dlmread(fname_neighborDistances,',',1,0);
    %neighborFeatures = neighborFeatures(:,2:(numNeighbors +1 ));
    
    fname_neighborIDs = fullfile(superpixel_input_dir,[im_splitStr{1},...
        '_se1_minNuc3_minStr5_minLum10_first50neigh_IDs']);
    neighborIDs = dlmread(fname_neighborIDs,',',1,0);
    
    fname_neighborRadii = fullfile(superpixel_input_dir,[im_splitStr{1},...
        '_se1_minNuc3_minStr5_minLum10_first50neigh_radii']);
    neighborRadii = dlmread(fname_neighborRadii,',',1,0);
    
    numNeighbors = size(neighborFeatures,2) - 1;
    numObjects = size(neighborFeatures,1);
    obj_radii = neighborRadii(:,1);
    
    % calculate the mean distances
    mean_neighborDists = mean(neighborDistances(:,2:(numNeighbors + 1)),2);
        
    objects_area = features_matrix(:,3);
    objects_radii = (objects_area./pi).^.5;
    figure; plot(obj_radii,mean_neighborDists,'b.');
    xlabel('Object Area'); ylabel('Mean neighbor distances');
    set(gca,'FontSize',16);
    figure; histogram(mean_neighborDists,'Normalization','probability','FaceColor',[0.8 0.8 0.8]);
    xlabel('Mean neighbor distance');
    set(gca,'FontSize',16);set(gcf,'color','white') 
    
    weights_neighbors = normpdf(neighborDistances(:,2:(numNeighbors + 1)),0,sigma)'; 
    weights_vectors = weights_neighbors(:);
    % convert the features matrix into pairs of pixels 
    F2 = neighborFeatures(:,2:(numNeighbors + 1))';
    F2 = F2(:);
    % object features
    object_features = repmat(neighborFeatures(:,1),[1 numNeighbors])';
    F1 = object_features(:);
    subplot(2,1,1);
    [histw, intervals] = histwc(features_matrix(:,2), features_matrix(:,3),50);
    bar(intervals, histw./sum(histw),'FaceColor',[0.8 .8 .8])
    ylim([0 0.2]);xlim([-pi pi]);set(gca,'FontSize',16);
    subplot(2,1,2);
    histogram(f_map(:),'Normalization','probability','FaceColor',[0.8 0.8 0.8],'NumBins',50);
    ylim([0 0.2]);xlim([-pi pi]);set(gca,'FontSize',16);

    % 2D weighted histogram
    cmap = colormap(linspecer(128));
    close(gcf);
    H = hist2w([F1 F2],weights_vectors./sum(weights_vectors),binCenters,binCenters,cmap);
    xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
    set(gcf,'color','white'); set(gca,'LooseInset',get(gca,'TightInset')); 

    figure; ndhist(F1(1:100:end,:),F2(1:100:end,:),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
    %figure; ndhist(F1,F2,'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
    xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
    set(gcf,'color','white'); set(gca,'LooseInset',get(gca,'TightInset'));

%     Y = zeros(numBins);
%     [bincount, indx_objs] = histc(neighborFeatures(:,1),binEdges);
%     
%     for indx = 1:numObjects
%         [bincount_neigh,indx_neigh] = histc(neighborFeatures(indx,2:end),binEdges);
%         weights_neighbor = normpdf(neighborDistances(indx,2:end),0,sigma);
%         for neigh = 1:numNeighbors
%             Y(indx_objs(indx),indx_neigh(neigh)) = Y(indx_objs(indx),indx_neigh(neigh)) ...
%                 + weights_neighbor(neigh);
%         end
%     end
%     c = [0.8 0.8 0.8];
%     figure;
%     b = bar3(Y); 
%     for k = 1:length(b)
%         zdata = b(k).ZData;
%         b(k).CData = zdata;
%         b(k).FaceColor = 'interp';
%     end
end







