% Demo file
% -------------------------------------------------------------------------
% HE tumor object segmentation 
% Luong Nguyen, lun5@pitt.edu
% Please email me if you find bugs, or have suggestions or questions

%% compile and check for error
addpath(genpath(pwd));
%compile;

%% Read in input file
% datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
% if ~ exist(datadir,'dir')
%     datadir = '/Users/lun5/Research/color_deconvolution/TissueImages/';
% end
datadir = 'test_images'; 
opts_input = setEnvironment_inputs;
% imname = 'gland3.tif';
% imname = 'gland1.tif';
% imname = 'gland3_snip.tif';
% imname = '101027.jpg'; % coral
% imname = '253027.jpg'; % zebra
% imname = '134067.jpg'; % leopard
% imname = 'fractal1.tif'; % fractal
% imname = 'tp10-611_22528_16384_2048_2048.tif';
imname = 'tp10-611gland7snip.tif';
% imname = 'tp10-867-1gland21.tif';
I = getImage(datadir, imname, opts_input);
% imshow(I); 
% rect = getrect; rect = round(rect)
% I = imcrop(I,rect);imshow(I);size(I)
%imwrite(I,fullfile('test_images','tp10-611gland7snip.tif'),'tif','Compression','none');
%% Calculate affinity matrix 
opts_affinity = setEnvironment_affinity;
[affinity_matrix, im_sizes] = calculateAffinity(I, opts_affinity);
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
