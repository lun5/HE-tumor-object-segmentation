% script to process the segmentation
tiles_dir = fullfile(pwd,'HEimages');
imname = '9uixINHtjjiS.tif';
%% result directory
splitStr = regexp(imname,'\.','split');
imresult_dir = fullfile(pwd,'results','HE_results',splitStr{1});

I = imread(fullfile(imresult_dir,'crop_image.tif'));
I = double(I);
% set environment
opts_clustering = setEnvironment_clustering;
opts_affinity = setEnvironment_affinity;
which_features = opts_affinity.features.which_features;
which_affinity = opts_affinity.affinityFunction;
methodresult_dir = fullfile(imresult_dir,[which_features{1} '_' which_affinity]);

% load E_oriented
E_oriented = load(fullfile(methodresult_dir,'E_oriented.mat'));
E_oriented = E_oriented.data;

if (~ispc)
        tic;thresh = 0.2;
        E_ucm = contours2ucm_crisp_boundaries(E_oriented,opts_affinity, opts_clustering);
        segmented_image = ucm2colorsegs(E_ucm,I,thresh);
        figure;
        subplot(121); imshow(uint8(I)); subplot(122); 
        imshow(uint8(segmented_image));toc
else
        segmented_image = [];
end