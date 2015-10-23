%% script to run Edison MShift Wrapper method
% Luong Nguyen 09/21/2015

github_dir = '/home/lun5/github/HE-tumor-object-segmentation';
mshift_dir = '/home/lun5/HEproject/evaluation_results/MeanShift';
im_dir = '/home/lun5/HEproject/data/Tiles_512';
result_dir = '/home/lun5/HEproject/evaluation_results/MeanShift';

% github_dir = '/Users/lun5/Research/github/HE-tumor-object-segmentation'; % mac
% mshift_dir = '/Users/lun5/Research/packages/edison_matlab_interface';
% im_dir = '/Users/lun5/Research/data/Tiles_512';
% result_dir = '/Users/lun5/Research/data/evaluation_results/MeanShift';

addpath(genpath(github_dir));
addpath(genpath(mshift_dir));

matfile_result_dir = fullfile(result_dir,'segmented_images');

if ~exist(result_dir,'dir')
    mkdir(result_dir)
end

if ~exist(matfile_result_dir,'dir')
    mkdir(matfile_result_dir);
end

if ~exist(fullfile(result_dir,'seg_im'),'dir')
    mkdir(fullfile(result_dir,'seg_im'));
end

if ~exist(fullfile(result_dir,'bdry_im'),'dir')
    mkdir(fullfile(result_dir,'bdry_im'));
end

im_list = dir(fullfile(im_dir,'*.tif'));
im_list = {im_list.name}';
% sigma, k, min
tmp = load('params_seism_mshift.mat'); params = tmp.params;
num_segs = size(params,1);
%run_times = cell(length(im_list),2);
parfor i = 1:length(im_list)
    im_name = im_list{i}(1:end-4);
    rgbim = imread(fullfile(im_dir,im_list{i}));
    segs = cell(num_segs,1);
    outFile = fullfile(matfile_result_dir,[im_name,'.mat']);
    if ~exist(outFile,'file')
        T = tic;
        for j = 1:size(params,1)
            % 'SpatialBandWidth'(sbw), 'RangeBandWidth'(rbw), 'MinimumRegionArea'(mra)
            sbw = params(j,1); rbw = params(j,2); mra = params(j,3);
            [~, labels, ~, ~, ~, ~] = edison_wrapper(rgbim, @RGB2Luv,...
                'SpatialBandWidth',sbw, 'RangeBandWidth', rbw,'MinimumRegionArea',mra);
            segs{j} = unit16(labels) + 1;
	    bdry_fname = fullfile(result_dir,'bdry_im',[im_name '_sbw' ...
                num2str(sbw) '_rbw' num2str(rbw) '_mra' num2str(mra) '.jpg']);
            seg_fname = fullfile(result_dir,'seg_im',[im_name '_sbw' ...
                num2str(sbw) '_rbw' num2str(rbw) '_mra' num2str(mra) '.jpg']);
            if ~ exist(bdry_fname,'file')
            	edge_map = seg2bdry(segs{j},'imageSize');
            	edge_map = imdilate(edge_map, strel('disk',1));
            	edge_map_im = rgbim.*uint8(repmat(~edge_map,[1 1 3]));
              	imwrite(edge_map_im,bdry_fname);
              	imwrite(label2rgb(segs{j}),seg_fname);
           end	
        end
        parsave(outFile,segs);
        t = toc(T);
        %run_times(i,:)= {im_name,t};
        fprintf('Done with image %s in %.2f s\n',im_name, t);        
    end
end

disp('Done');

