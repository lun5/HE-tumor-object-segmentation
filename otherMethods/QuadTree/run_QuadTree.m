%% script to run QuadTree method
% Luong Nguyen 09/21/2015
% Code provided by Jordi Pont-Tuset

% github_dir = '/home/lun5/github/HE-tumor-object-segmentation';
% mshift_dir = '/home/lun5/HEproject/evaluation_results/QuadTree';
% im_dir = '/home/lun5/HEproject/data/Tiles_512';
% result_dir = '/home/lun5/HEproject/evaluation_results/QuadTree';

github_dir = '/Users/lun5/Research/github/HE-tumor-object-segmentation'; % mac
QuadTree_dir = '/Users/lun5/Research/github/seism/src/baseline';
im_dir = '/Users/lun5/Research/data/Tiles_512';
result_dir = '/Users/lun5/Research/data/evaluation_results/QuadTree';

addpath(genpath(github_dir));
addpath(genpath(QuadTree_dir));

matfile_result_dir = fullfile(result_dir,'segmented_images');

if ~exist(result_dir,'dir')
    mkdir(result_dir)
end

if ~exist(matfile_result_dir,'dir')
    mkdir(matfile_result_dir);
end

im_list = dir(fullfile(im_dir,'*.tif'));
im_list = {im_list.name}';
% sigma, k, min
tmp = load('../params_seism_QuadTree.mat'); params = tmp.params;
num_segs = size(params,1);
%run_times = cell(length(im_list),2);
sx = 512; sy = 512; % image size. One segmentations for all images
grid_ucm = create_quad_tree(sx,sy);
segs = cell(num_segs,1);
parfor i = 1:num_segs
    nseg = params(i);
    segs{i} = griducm2seg(grid_ucm,nseg);
end

parfor i = 1:length(im_list)
    im_name = im_list{i}(1:end-4);
    outFile = fullfile(matfile_result_dir,[im_name,'.mat']);
    if ~exist(outFile,'file')
        parsave(outFile,segs);
    end
end

disp('Done');

