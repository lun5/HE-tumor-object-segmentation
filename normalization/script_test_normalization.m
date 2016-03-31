%% test script for normalization
%sourcedir = 'Z:\';
%tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%tiles_dir = 'Z:\Tiles_512';
tiles_dir = 'C:\Users\luong_nguyen\Box Sync\ADH\ImageColorNormalization';
%% test for the ADH code
target_im_name = 'source1.tif';
source_im_name = 'target2.tif';

source_im = imread(fullfile(tiles_dir,source_im_name));
%target_im = imread(fullfile(tiles_dir,target_im_name));
target_im = imread('Z:\Tiles_512\2ale5ngryfnpo.tif');
%target_im = imread('Z:\Tiles_512\ocmmhhrtzz5.tif');
rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
opts = setEnvironment_affinity;
opts.matchMethod = 'moments';

mask_source = ~isolateWhite(source_im);
%[ source_eq_image, ~,~,~ ] = ...
%    oppColNormalization( source_im, target_im, rotation_matrix, opts);
tic;
[ source_eq_image] = col_normalization( source_im, mask_source);
toc
imwrite(source_eq_image,fullfile(tiles_dir,['new_eq_' source_im_name]));
%%
%out_dir =  fullfile(sourcedir,'HEproject','normalization_cl_isolateRW');
% out_dir =  fullfile(sourcedir,'HEproject','normalization_512_target95f7');
% %mV_dir = 'Z:\mixture_von_mises\isolate_white_red_nonfreeze';
% %mV_dir = 'Z:\mixture_von_mises\isolate_white_red_nonfreeze_restrict_purple_spread';
% mV_dir = fullfile(sourcedir,'mixture_von_mises','isolate_white_red_purple_512');
% 
% if ~exist(out_dir,'dir')
%     mkdir(out_dir);
% end
% target_im_name = '95f7k8LoesYevi.tif';
% %source_im_name = 'P3msE3FrHJz.tif';
% %source_im_name = 'pBPHL1xUjdvYx.tif';
% %source_im_name = '6NfelbSiNWJxYst.tif';
% %target_im_name = 'oCMMhhrTzZ5.tif';
% fileNames = dir(fullfile(tiles_dir,'*.tif'));
% imagepaths = {fileNames.name}';
% numImages = length(imagepaths);% 232
% 
% %source_im = imread(fullfile(tiles_dir,source_im_name));
% %target_im = imread(fullfile(tiles_dir,target_im_name));
% %figure; imshow(source_im);
% %figure; imshow(target_im);
% 
% rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
% opts = setEnvironment_affinity;
% opts.matchMethod = 'moments'; 
% parfor i = 1:numImages
%     source_im_name = imagepaths{i};
%     if exist(fullfile(out_dir,source_im_name),'file') || strcmp(source_im_name,target_im_name)
%         continue;
%     end
%     t = tic;
%     [ source_eq_image, f_maps_source, f_maps_target,f_maps_source_normalized ] ...
%         = opp_col_normalization( source_im_name(1:end-4), target_im_name(1:end-4),...
%         rotation_matrix, tiles_dir, mV_dir, out_dir, opts);
%     T = toc(t); fprintf('done with image %s in %.2f seconds\n',source_im_name(1:end-4),T);
%     close all;
% end
% disp('Done');


%% samples of pair of pixels
% how do the histogram look different 
% Nsamples = 10000;
% [source_sample,p1,p2] = sampleF(f_maps_source{1},Nsamples,opts);
% [target_sample,p1,p2] = sampleF(f_maps_target{1},Nsamples,opts);
% 
% Fsym_source = [source_sample; source_sample(:,2) source_sample(:,1)];
% figure; ndhist(Fsym_source(:,1),Fsym_source(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
% xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
% set(gcf,'color','white') 
% 
% Fsym_target = [target_sample; target_sample(:,2) target_sample(:,1)];
% figure; ndhist(Fsym_target(:,1),Fsym_target(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
% xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
% set(gcf,'color','white') 
% 
% [source_eq_sample,p1,p2] = sampleF(f_maps_source_normalized{1},Nsamples,opts);
% Fsym_source_eq = [source_eq_sample; source_eq_sample(:,2) source_eq_sample(:,1)];
% figure; ndhist(Fsym_source_eq(:,1),Fsym_source_eq(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
% xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
% set(gcf,'color','white')