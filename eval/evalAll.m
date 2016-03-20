%% function [] = evalAll(IMG_DIR,GT_DIR,RESULTS_DIR,type)
%
% INPUTS
%  IMG_DIR     - path to images (assumes images are jpegs)
%  GT_DIR      - path to ground truth boundary and segmentation maps
%  RESULTS_DIR - where to write out the results
%
% OUTPUTS
%  none, will write results to RESULTS_DIR and plot the metrics
% 
% See evaluation_script.m for an example of how to call this function
%
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2015 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------


function [] = evalAll(IMG_DIR,GT_DIR,RESULTS_DIR, opts_affinity, opts_clustering)
    
    %% read images
    IMG_EXT = '.tif';
    img_list = dirrec(IMG_DIR,IMG_EXT);
    fprintf('Number of images is %d',length(img_list));
    %% compute boundaries for images
    if (~exist(RESULTS_DIR,'dir'))
        mkdir(RESULTS_DIR);
    end
    if (~exist(fullfile(RESULTS_DIR,'E_oriented'),'dir'))
        mkdir(fullfile(RESULTS_DIR,'E_oriented'));
    end
    if (~exist(fullfile(RESULTS_DIR,'ucm2'),'dir'))
        mkdir(fullfile(RESULTS_DIR,'ucm2'));
    end
    if (~exist(fullfile(RESULTS_DIR,'segmented_images'),'dir'))
        mkdir(fullfile(RESULTS_DIR,'segmented_images'));
    end
    if (~exist(fullfile(RESULTS_DIR,'edgemap'),'dir'))
        mkdir(fullfile(RESULTS_DIR,'edgemap'));
    end
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_all232');
    EV_DIR = fullfile(RESULTS_DIR,'ev_txt_fine');
    if (~exist(EV_DIR,'dir'))
        mkdir(EV_DIR);
    end    
    %opts_affinity = setEnvironment_affinity;
    
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});im_name = lower(im_name);
        if (~exist(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),'file'))
            fprintf('\n\nCalculate E oriented %s...',im_name); tic;
            I = imread(img_list{i});
            I = double(I);
	    mult = 1; I = I(1:mult:end,1:mult:end,:);
            [Ws,Ws_each_feature_set, im_sizes] = getW(I,opts_affinity);
            [~, ~,~, E_oriented] = graphSegmentation(Ws,Ws_each_feature_set{1},im_sizes,I,opts_clustering);
            E_oriented = imresize(E_oriented{1},size(I(:,:,1)));
            E = max(E_oriented,[],3);
            imwrite(mat2gray(1-E),fullfile(RESULTS_DIR,'edgemap',[im_name,'_edgemap.tif']),'Resolution',300); 
            parsave(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),E_oriented);
            t = toc; fprintf('done: %1.2f sec\n', t);
        end      
    end
   
    %% normalize output scale
    E_orienteds = cell(1,length(img_list));
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});im_name = lower(im_name);
        tmp = load(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']));
        E_orienteds{i} = tmp.data;
    end
    %
    max_val = max(cellfun(@max_all,E_orienteds));
    parfor i=1:length(img_list)
        E_orienteds{i} = E_orienteds{i}/max_val;
    end
    %% run UCM on boundary maps
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});im_name = lower(im_name);
        if (~exist(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']),'file'))
            fprintf('\n\nCalculate UCM %s...',im_name); T = tic;
            ucm2 = contours2ucm_crisp_boundaries(mat2gray(E_orienteds{i}),'doubleSize');
            parsave(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']),ucm2);
            t = toc(T); fprintf('done: %1.2f sec\n', t);
        end
%         if ~exist(fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif']),'file')
%             fprintf('\n\nCalculate segmented image %s...',im_name); T = tic;
%             tmp = load(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']));
%             ucm2 = tmp.data;
%             I = imread(img_list{i});I = double(I(1:4:end,1:4:end,:));
%             segmented_image = ucm2colorsegs(ucm2,I,0.2);
%             imwrite(uint8(segmented_image),fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif']),'Resolution',300); 
%             t = toc(T); fprintf('done: %1.2f sec\n', t);
%         end
    end
    
    %% eval using BSR metrics
    %allBench_custom(IMG_DIR,GT_DIR,fullfile(RESULTS_DIR,'ucm2'),EV_DIR);
    %eval_Fop(IMG_DIR, GT_DIR, fullfile(RESULTS_DIR,'ucm2'),EV_DIR); 
    %plot_eval(EV_DIR);
end
