%% function [] = evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR)
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


function [] = evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR)
    
    %% read images
    %IMG_EXT = '.tif';
    %img_list = dirrec(IMG_DIR,IMG_EXT);
    %% compute boundaries for images
    if (~exist(RESULTS_DIR,'dir'))
        mkdir(RESULTS_DIR);
    end
    
    SEG_DIR = fullfile(RESULTS_DIR,'segmented_images');
    img_list = dirrec(SEG_DIR,'.mat');
    fprintf('Number of files is %d\n',length(img_list));
    EV_DIR = fullfile(RESULTS_DIR,'ev_txt');
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_all232');
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_reannotated_Oct14');
    if (~exist(EV_DIR,'dir'))
        mkdir(EV_DIR);
    end
   
    %[~,im_name,~] = fileparts(img_list{1});    
    %tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    %segs = tmp.data; 
    %nSegments = length(segs); % segments 2:2:200
    nSegments = 1;
    %% eval using BSR metrics
    allBench_custom(IMG_DIR,GT_DIR,SEG_DIR,EV_DIR,nSegments);
    eval_Fop(IMG_DIR, GT_DIR, SEG_DIR, EV_DIR,nSegments);
    plot_eval(EV_DIR);


