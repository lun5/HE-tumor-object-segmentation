% Luong Nguyen 2/29/2016
%% connected component test
GIT_DIR = '/home/lun5/github/tissue-component-classification';
addpath(genpath(GIT_DIR));
%IMG_DIR = 'Z:\ADH_Jeff\CVPR_images\Breast WSI object lists';
%IMG_DIR ='/home/lun5/CVPR_images/breast_wsi';
%IMG_DIR = '/home/lun5/HEproject/CVPR_images';
%/home/lun5/CVPR_images/breast_wsi/adjDela
%IMG_DIR = 'C:\Users\luong_nguyen\Box Sync\ADH\CVPRImages Superpixels\Breast WSI object lists';
%imlist = dir(fullfile(IMG_DIR,'images','*.jpg'));
%imlist = dir(fullfile(IMG_DIR,'images','*.tif'));
%IMG_DIR = 'Z:\ADH_Jeff\CVPR_images\';
IMG_DIR = 'Z:\HEproject\object_proposals';
%imlist = dir(fullfile(IMG_DIR,'adjDela','*_se1_*'));
imlist = dir(fullfile(IMG_DIR,'AdjDela_updated','*_se1_*'));
imlist = {imlist.name}';
%imlist = {'tp10-867-1','tp10-876-1'}
%imlist = {'1050508_ImMod276_stX22001_stY24001.jpg'};
num_images = length(imlist);
%param_string = '_se1_minNuc3_minStr5_minLum9';
%param_string = '_se1_minNuc5_minStr5_minLum9';
param_string = '_se1_minNuc3_minStr5_minLum5';
%param_string ='_se1_minNuc3_minStr5_minLum9'
num_neighbors = 15;
num_comps = 20; 
plot_flag = 3;
top_centers = cell(num_images,1);
top_radii = cell(num_images,1);
im_size = [2048 2048];
% output directory
%outdir = fullfile(IMG_DIR, ['mutual_knn_cc_' num2str(num_neighbors)]);
outdir = fullfile(IMG_DIR,'updated_cca_voronoi');
if ~exist(outdir, 'dir');
    mkdir(outdir);
end

if ~ exist(fullfile(outdir,'segmented_images'),'dir')
    mkdir(fullfile(outdir,'segmented_images'));
end

if ~ exist(fullfile(outdir,'bdry_im'),'dir')
    mkdir(fullfile(outdir,'bdry_im'));
end

for i = 1:num_images
    T = tic;
    imname = imlist{i}(1:end-36);
    %imname = '9uixINHtjjiS';
    %imname = '2ale5ngryfnpo';
    %imname = '1yDJEq9l1cuUNk7';
    %imname = 'HlN7HDXf2aFQ';
    %imname = 'YXDtQ9VZRKUE1';
    %imname = 'fFwTGXYlhYNa';
    imname = 'lszomRlGsC5na4Q';
%     if exist(fullfile(outdir,'segmented_images',[imname '.mat']),'file')
%         continue;
%     end
    segs = object_proposal_all_types( IMG_DIR, imname, param_string, ...
        num_neighbors, num_comps, plot_flag ); 
    runtime = toc(T);
    fprintf('Done with %s in %.2f seconds.\n',imname,runtime);
    %close all;
end

% imname = 'tp09-1003';
% num_comp = 100; num_neighbors = 15;
% [top_centers, top_radii] = top_connected_comp( IMG_DIR, imname, param_string, ...
%     obj_type, num_neighbors, num_comps, plot_flag ); 
disp('Done');

% I = double(imread(fullfile(IMG_DIR, [imname '.jpg'])));
% mask = dlmread(fullfile(IMG_DIR,[imname '_k3']),',',1,0);
% mask_purple = uint8(mask == 1);
% 
% circle_map = dlmread(fullfile(IMG_DIR,[imname '_se1_minNuc3_minStr5_minLum9_circle_map']),',',1,0);
% adj_map = dlmread(fullfile(IMG_DIR,[imname '_se1_minNuc3_minStr5_minLum9_adjDela']),',',0,0);
% obj_types = adj_map(:,3);
% obj_coords = adj_map(:,[5 4]);
% obj_radii = sqrt(adj_map(:,2)./pi);
% 
% obj_type = 1; % nuclei 1, stroma 2, lumen 3
% indx_type = obj_types == obj_type;
% centers = obj_coords(indx_type,:);
% radii = obj_radii(indx_type);
% obj_map = mat2gray(circle_map == obj_type);
% 
% figure; imshow(obj_map);
% hold on; viscircles(centers, radii);
% 
% %% Now test with LOG
% %% dilated mask
% se = strel('disk',2,4);
% dilate_mask = imdilate(mask_purple,se);
% masked_I = I./255.*double(repmat(mask_purple,[1 1 3]));
% masked_im = mean(masked_I,3);
% masked_im(mask ~= 1) = 1; figure; imshow(masked_im);
% % LOG parameters
% octave = 2; intervals = 2;
% contrast_threshold = 0.01;
% curvature_threshold = 10.0;
% initial_sigma = 4;
% max_choice = 1;
% area_threshold = 0.1;
% [centers, radii] = scale_orient_selection(masked_im(1:500,1:500), octave, intervals,...
%    dilate_mask(1:500,1:500), contrast_threshold, curvature_threshold,initial_sigma,max_choice, area_threshold, 1);
% 
% % euclidian distance between centers 
% D = pdist(centers,'euclidean');
% %figure; histogram(D);
% % threshold to filter out pairs of pixels that are not too close by
% dist_thres = 60; % pixels for now
% indx = D > dist_thres;
% D( indx ) = 0;
% sigma = std(D); % standard deviation of the distance
% similarities = exp(- D.^2./sigma);
% similarities(indx) = 0;
% A = squareform(similarities);
% %G = graph(uint8(A>0));
% %plot(G);
% 
% [S,C] = graphconncomp(sparse(A));
% num_comps = 10;
% colors = distinguishable_colors(num_comps);
% figure; imshow(I);
% for i = 1:num_comps
%     indx_cl = find(C == i+20);
%     viscircles(centers(indx_cl,:), radii(indx_cl),'EdgeColor',colors(i,:));
% end

% LOG is not doing so well with connected components. 

%% connected component for stroma and lumen
