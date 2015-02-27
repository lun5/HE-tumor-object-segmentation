% Demo file
% -------------------------------------------------------------------------
% HE tumor object segmentation 
% Luong Nguyen, lun5@pitt.edu
% Please email me if you find bugs, or have suggestions or questions

%% compile and check for error
%addpath(genpath(pwd));
%compile;

%% Read in input file
% datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
% if ~ exist(datadir,'dir')
%     datadir = '/Users/lun5/Research/color_deconvolution/TissueImages/';
% end
%sourcedir = 'Z:\';
%tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
tiles_dir = '/Users/lun5/Box Sync/TilesForLabeling_tiff_renamed';
opts_input = setEnvironment_inputs;

%I = getImage(datadir, imname, opts_input);
%I = imread(fullfile(tiles_dir, 'fFwTGXYlhYNa.tif'));
%I = imread(fullfile(tiles_dir, '5aOQp0sbfXSWMZ.tif'));
I = imread(fullfile(tiles_dir, '9uixINHtjjiS.tif'));
%I = imread(fullfile(tiles_dir, '46vr2niG4yne5Lx.tif'));
%I = imread(fullfile(tiles_dir, 'BQC7vv3HUhCe.tif'));
%I = imread(fullfile(tiles_dir, 'P3msE3FrHJz.tif'));
%I = imread(fullfile(tiles_dir, 'EMnOxgxqoMGzn1.tif'));
%I = imread(fullfile(tiles_dir, 'ApAaL7fc2paYi.tif'));
%I = imread(fullfile(tiles_dir, 'IlGwtTFXmQ.tif'));
%I = imread(fullfile(tiles_dir, 'dJUtEn6DHnfd.tif'));
%I = imread(fullfile(pwd,'test_images','253027.jpg'));
%I = imread(fullfile(pwd,'fractal','fracTest1.pgm'));
%I = imread(fullfile(tiles_dir, 'tp09-96-2_10240_28672_2048_2048.tif'));


% [I, segIm] = rbfFracImageNew([],[],[],[50 50]); 
close all; 
clear E_oriented opts_affinity opts_clustering segmented_image affinity_matrix
% rect = getrect; rect = round(rect);
% I = imcrop(I,rect);imshow(I);size(I)
%imwrite(I,fullfile('test_images','tp10-611gland7snip.tif'),'tif','Compression','none');
%% Calculate affinity matrix 
I_downsample = imresize(I,1/4);figure;imshow(I_downsample);
opts_affinity = setEnvironment_affinity;
[affinity_matrix, im_sizes] = calculateAffinity(I_downsample, opts_affinity);
% need to find a place to put the rotation matrix in somewhere

%% Graph-based clustering based on 
% this depends on whether the outputs are segmentation or detecting edges
opts_clustering = setEnvironment_clustering;
[segmented_image, E_oriented] = graphSegmentation(affinity_matrix,im_sizes,I,opts_clustering,opts_affinity);

% %% Display edges and segmentation
% % get gland dat afor tp10-611
% svs_fname = 'tp10-611';
% datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
% if ~ exist(datadir,'dir')
%     datadir = '/Users/lun5/Research/color_deconvolution/TissueImages/';
% end
% resultdir = fullfile(pwd,'test_images');
% wsi_get_objects(svs_fname, datadir, resultdir)
