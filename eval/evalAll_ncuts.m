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
    
    nSegments = 99;
    parfor i = 1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        fprintf('\n\nCalculate Ncuts Segmentation %s...',im_name);T = tic;
        outFile = fullfile(RESULTS_DIR,'segmented_images',[im_name '.mat']);
        I = double(imread(img_list{i}));
        I = mean(I,3); [nr,nc] = size(I);
        seg = cell(nSegments,1);
        %mult = 4; I = double(I(1:mult:end,1:mult:end,:));    
        [W,~] = ICgraph(I);
        [NcutEigenvectors,~] = ncut(W,nSegments + 1);
        for j= 2:nSegments+1
            %[NcutDiscrete,~,~] = ncutW(W,j);
            [NcutDiscrete,~] = discretisation(NcutEigenvectors(:,1:j));
            SegLabel = zeros(nr,nc);
            for k=1:size(NcutDiscrete,2),
                SegLabel = SegLabel + k*reshape(NcutDiscrete(:,k),nr,nc);
            end            
            seg{j-1} = SegLabel;            
        end
        t = toc(T); fprintf('done: %1.2f sec\n', t);
        parsave(outFile,seg);
    end
    
    %% eval using BSR metrics
    SEG_DIR = fullfile(RESULTS_DIR,'segmented_images');
    %allBench_custom(IMG_DIR,GT_DIR,SEG_DIR,fullfile(RESULTS_DIR,'ev_txt'),nSegments);
    %plot_eval(fullfile(RESULTS_DIR,'ev_txt'));

