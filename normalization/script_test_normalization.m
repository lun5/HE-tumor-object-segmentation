%% test script for normalization
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
out_dir =  fullfile(sourcedir,'HEproject','normalization_cl');
mV_dir = 'Z:\mixture_von_mises\isolate_white_nonfreeze';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
%source_im_name = '95f7k8LoesYevi.tif';
source_im_name = 'P3msE3FrHJz.tif';
%source_im_name = 'pBPHL1xUjdvYx.tif';
%source_im_name = '6NfelbSiNWJxYst.tif';
target_im_name = 'oCMMhhrTzZ5.tif';
fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths);% 232

%source_im = imread(fullfile(tiles_dir,source_im_name));
%target_im = imread(fullfile(tiles_dir,target_im_name));
%figure; imshow(source_im);
%figure; imshow(target_im);

rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
opts = setEnvironment_affinity;
opts.matchMethod = 'moments'; 
for i = 1:numImages
    source_im_name = imagepaths{i};
    tic;
    [ source_eq_image, f_maps_source, f_maps_target,f_maps_source_normalized ] ...
        = opp_col_normalization( source_im_name(1:end-4), target_im_name(1:end-4),...
        rotation_matrix, tiles_dir, mV_dir, out_dir, opts);
    T = toc; fprintf('done with image %s in %.2f\n',source_im_name(1:end-4),T);
    close all;
end
%% samples of pair of pixels
% how do the histogram look different 
Nsamples = 10000;
[source_sample,p1,p2] = sampleF(f_maps_source{1},Nsamples,opts);
[target_sample,p1,p2] = sampleF(f_maps_target{1},Nsamples,opts);

Fsym_source = [source_sample; source_sample(:,2) source_sample(:,1)];
figure; ndhist(Fsym_source(:,1),Fsym_source(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 

Fsym_target = [target_sample; target_sample(:,2) target_sample(:,1)];
figure; ndhist(Fsym_target(:,1),Fsym_target(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 

[source_eq_sample,p1,p2] = sampleF(f_maps_source_normalized{1},Nsamples,opts);
Fsym_source_eq = [source_eq_sample; source_eq_sample(:,2) source_eq_sample(:,1)];
figure; ndhist(Fsym_source_eq(:,1),Fsym_source_eq(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white')