% Demo file
% -------------------------------------------------------------------------
% HE tumor object segmentation 
% Luong Nguyen, lun5@pitt.edu
% Please email me if you find bugs, or have suggestions or questions

%% compile and check for error
addpath(genpath(pwd));
compile;

%% Read in input file
datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
if ~ exist(datadir,'dir')
    datadir = '/Users/lun5/Research/color_deconvolution/TissueImages/';
end

opts_input = setEnvironment_inputs();
imname = 'tp10-867-1_47104_22528_2048_2048.tif';
I = getImage(datadir, imname, opts_input);

%% Calculate affinity matrix 
opts_affinity = setEnvironment_affinity();
affinity_matrix = calculateAffinity(I, opts);
% need to find a place to put the rotation matrix in somewhere

%% Graph-based clustering based on 
% this depends on whether the outputs are segmentation or detecting edges
opts_clustering = setEnvironment_clustering();
[segmented_image, numComponents] = graphSegmentation(I, affinity_matrix, opts_clustering);

%% Display edges and segmentation


%% Segment image
% builds an Ultrametric Contour Map from the detected boundaries (E_oriented)
% then segments image based on this map
%
% this part of the code is only supported on Mac and Linux

if (~ispc)
    
    thresh = 0.1; % larger values give fewer segments
    E_ucm = contours2ucm_crisp_boundaries(E_oriented,type);
    S = ucm2colorsegs(E_ucm,I,thresh);

    close all; subplot(121); imshow(I); subplot(122); imshow(S);
end
