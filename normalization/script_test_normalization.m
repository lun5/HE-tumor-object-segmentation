%% test script for normalization
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%source_im_name = '95f7k8LoesYevi.tif';
source_im_name = 'P3msE3FrHJz.tif';
%source_im_name = 'pBPHL1xUjdvYx.tif';
%source_im_name = '6NfelbSiNWJxYst.tif';
target_im_name = 'oCMMhhrTzZ5.tif';

source_im = imread(fullfile(tiles_dir,source_im_name));
target_im = imread(fullfile(tiles_dir,target_im_name));
figure; imshow(source_im);
figure; imshow(target_im);

rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
opts = setEnvironment_affinity;
opts.matchMethod = 'moments'; 
tic;
[ source_eq_image, f_maps_source, f_maps_target,f_maps_source_normalized ] ...
     = oppColNormalization( source_im, target_im,rotation_matrix, opts);
toc
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