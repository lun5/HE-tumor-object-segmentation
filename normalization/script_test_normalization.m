%% test script for normalization
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%source_im_name = '95f7k8LoesYevi.tif';
%source_im_name = 'P3msE3FrHJz.tif';
%source_im_name = 'pBPHL1xUjdvYx.tif';
source_im_name = '6NfelbSiNWJxYst.tif';
target_im_name = 'oCMMhhrTzZ5.tif';

source_im = imread(fullfile(tiles_dir,source_im_name));
target_im = imread(fullfile(tiles_dir,target_im_name));
figure; imshow(source_im);
figure; imshow(target_im);

rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
opts = setEnvironment_affinity;
opts.matchMethod = 'moments'; 
tic;
[ source_eq_image ] = oppColNormalization( source_im, target_im,rotation_matrix, opts);
toc