%% function [] = evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR)
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


function [] = evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR,ev_mode)
    
    %% read images
    %IMG_EXT = '.tif';
    %img_list = dirrec(IMG_DIR,IMG_EXT);
    
    GT_DIR_mode = fullfile(GT_DIR,ev_mode);
    IMG_EXT = '.mat';
    img_list = dirrec(GT_DIR_mode,IMG_EXT);
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
    if (~exist(fullfile(RESULTS_DIR,'montageImageSegment'),'dir'))
        mkdir(fullfile(RESULTS_DIR,'montageImageSegment'));
    end

    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_all232');
    EV_DIR = fullfile(RESULTS_DIR,['ev_txt_' ev_mode]);
    if (~exist(EV_DIR,'dir'))
        mkdir(EV_DIR); 
    end
    % note that bsr only take image name   
   
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});%im_name = lower(im_name); 
        outFile = fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']);
        fprintf('\n\nCalculate E oriented %s...',im_name); T = tic;
        if (~exist(outFile,'file'))
            %E_oriented = globalPb(img_list{i}, outFile);
            E_oriented = proxy_globalPb(img_list{i}, outFile);
            parsave(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),E_oriented);
        end
	if (~exist(fullfile(RESULTS_DIR,'edgemap',[im_name,'_edgemap.tif']),'file')) 
	    tmp = load(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']));
            E_oriented = tmp.data;
	    E = max(E_oriented,[],3);
            imwrite(mat2gray(1-E),fullfile(RESULTS_DIR,'edgemap',[im_name,'_edgemap.tif']),'Resolution',300); 
            t = toc(T); fprintf('done: %1.2f sec\n', t);         
        end
    end
    
    %% normalize output scale
    E_orienteds = cell(1,length(img_list));
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        tmp = load(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']));
        E_orienteds{i} = tmp.data;
    end
    
    max_val = max(cellfun(@max_all,E_orienteds));
    parfor i=1:length(img_list)
        E_orienteds{i} = E_orienteds{i}/max_val;
    end
    
    %% run UCM on boundary maps
    UCM_DIR = fullfile(RESULTS_DIR,'ucm2');
    %{
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        if (~exist(fullfile(UCM_DIR,[im_name '.mat']),'file'))
            fprintf('\n\nCalculate UCM %s...',im_name); T = tic;
            ucm2 = contours2ucm(mat2gray(E_orienteds{i}), 'doubleSize');            
    %        ucm2 = proxy_contours2ucm(mat2gray(E_orienteds{i}),'doubleSize');
            parsave(fullfile(UCM_DIR,[im_name '.mat']),ucm2);
            t = toc(T); fprintf('done: %1.2f sec\n', t);
        end        
    end
    %}    
    %% eval using BSR metrics
    allBench_custom(IMG_DIR,GT_DIR_mode,UCM_DIR,EV_DIR);
    eval_Fop(IMG_DIR, GT_DIR_mode, UCM_DIR, EV_DIR);
    plot_eval(EV_DIR);

