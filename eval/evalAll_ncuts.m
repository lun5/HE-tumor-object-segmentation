%% function [] = evalAll_ncuts(IMG_DIR,GT_DIR,RESULTS_DIR)
% INPUTS
%  IMG_DIR     - path to images (assumes images are jpegs)
%  GT_DIR      - path to ground truth boundary and segmentation maps
%  RESULTS_DIR - where to write out the results
%
% OUTPUTS
%  none, will write results to RESULTS_DIR and plot the metrics
% 
% See evaluation_script.m for an example of how to call this function
% function to evaluate performance of BSR
% based on bsr code + Isola's code
% Luong Nguyen 8/6/2015
% lun5@pitt.edu


function [] = evalAll_ncuts(IMG_DIR,GT_DIR,RESULTS_DIR)
    
    %% read images
    IMG_EXT = '.tif';
    img_list = dirrec(IMG_DIR,IMG_EXT);

    %% compute boundaries for images
    if (~exist(RESULTS_DIR,'dir'))
        mkdir(RESULTS_DIR);
    end
    if (~exist(fullfile(RESULTS_DIR,'segmented_images'),'dir')) %seg label
        mkdir(fullfile(RESULTS_DIR,'segmented_images'));
    end
    if (~exist(fullfile(RESULTS_DIR,'ev_txt'),'dir'))
        mkdir(fullfile(RESULTS_DIR,'ev_txt'));
    end    
    
    nSegments = 100;
    for i = 1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        fprintf('\n\nCalculate E oriented %s...',im_name); tic;
        outFile = fullfile(RESULTS_DIR,'segmented_images',[im_name '_SegLabel.mat']);
        I = imread(img_list{i});
        I = mean(I,3);
        seg = cell(nSegments,1);
        mult = 4; I = double(I(1:mult:end,1:mult:end,:));      
        parfor j=1:nSegments
            [SegLabel,~,~,~,~,~]= NcutImage(I,j);
            seg{j} = SegLabel;
        end
        parsave(outFile,seg);
    end
    
    %% eval using BSR metrics
    SEG_DIR = fullfile(RESULTS_DIR,'segmented_images');
    allBench_custom(IMG_DIR,GT_DIR,SEG_DIR,fullfile(RESULTS_DIR,'ev_txt'));
    plot_eval(fullfile(RESULTS_DIR,'ev_txt'));

