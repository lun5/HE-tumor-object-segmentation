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
% Luong modified on 7/10/15 to combine ucm

function [] = evalAll_new(IMG_DIR,GT_DIR,RESULTS_DIR, opts_affinity)
    
    %% read images
    IMG_EXT = '.tif';
    img_list = dirrec(IMG_DIR,IMG_EXT);
    
    %% compute boundaries for images
    if (~exist(RESULTS_DIR,'dir'))
        mkdir(RESULTS_DIR);
    end
    if (~exist(fullfile(RESULTS_DIR,'E_oriented'),'dir'))
        mkdir(fullfile(RESULTS_DIR,'E_oriented'));
    end
    UCM_DIR = fullfile(RESULTS_DIR,'ucm2');
    if (~exist(UCM_DIR,'dir'))
        mkdir(UCM_DIR);
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
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt');
    EV_DIR = fullfile(RESULTS_DIR,'ev_txt_all232');
    if (~exist(EV_DIR,'dir'))
        mkdir(EV_DIR);
    end    
    
    %opts_affinity = setEnvironment_affinity;
    opts_clustering = setEnvironment_clustering;
    opts_clustering.display_progress = false;
    opts_clustering.calculate_segments = false;
    opts_clustering.plot_results = false;  
    mult = 1; % subsampling the image
       
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        if (~exist(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),'file'))
            fprintf('\n\nCalculate E oriented %s...',im_name); tic;
            I = imread(img_list{i});
            dz_im = I(1:mult:end,1:mult:end,:);
            I = double(dz_im);
            [Ws,Ws_each_feature_set, im_sizes] = getW(I,opts_affinity);
            [~,~, E, E_oriented] = graphSegmentation(Ws,Ws_each_feature_set{1},im_sizes,I,opts_clustering);
            imwrite(mat2gray(1-E),fullfile(RESULTS_DIR,'edgemap',[im_name,'_edgemap.tif']),'Resolution',300); 
            parsave(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),E_oriented);
            t = toc; fprintf('done: %1.2f sec\n', t);         
        end
    end
    
    %% normalize output scale
    E_orienteds = cell(1,length(img_list));
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        tmp = load(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']));
        E_orienteds{i} = tmp.data;
        max_vals{i} = max(cellfun(@max_all,E_orienteds{i}));
    end
    %
    max_val = max(cellfun(@max_all,max_vals));
    parfor i=1:length(img_list)
        for j = 1:length(E_orienteds{i})
            tmp = E_orienteds{i}{j};
            E_orienteds{i}{j} = tmp./max_val;
        end
        %E_orienteds{i} = E_orienteds{i}/max_val;
    end    
    
    %% run UCM on boundary maps
    weights = [10 2 1]'; % weights to combine hue, brightness, saturation
    UCM_DIR = fullfile(RESULTS_DIR,'ucm2',['weights_', num2str(weights','%d_%d_%d')]);
    if (~exist(UCM_DIR,'dir'))
        mkdir(UCM_DIR)
    end
    
    parfor i =1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        if (~exist(fullfile(UCM_DIR,[im_name '.mat']),'file'))
            fprintf('\n\nCalculate UCM %s...',im_name); T1 = tic;
            E_oriented = E_orienteds{i}; num_channels = length(E_oriented);
            all_ucms = cell(1,num_channels);
            for channel = 1:num_channels
                all_ucms{channel} = contours2ucm_crisp_boundaries(mat2gray(E_oriented{channel}),'doubleSize');
            end
            all_ucms = cat(3,all_ucms{:});
            W = repmat(weights./sum(weights),1, size(all_ucms,1), size(all_ucms,2));
            W = permute(W,[2 3 1]);
            ucm2 = sum(all_ucms.*W,3);
            parsave(fullfile(UCM_DIR,[im_name '.mat']),ucm2);
            t1 = toc(T1); fprintf('done: %1.2f sec\n', t1);
        end
%         if ~exist(fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif']),'file')
%             fprintf('\n\nCalculate segmented image %s...',im_name); T2=tic;
%             tmp = load(fullfile(UCM_DIR,[im_name '.mat']));
%             ucm2 = tmp.data;
%             I = imread(img_list{i});I = double(I(1:4:end,1:4:end,:));
%             [segmented_image,~] = ucm2colorsegs(ucm2,I,0.2);
%             imwrite(uint8(segmented_image),fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif']),'Resolution',300); 
%             t2 = toc(T2); fprintf('done: %1.2f sec\n', t2);
%         end
    end
    E_orienteds = [];
    %% eval using BSR metrics
    allBench_custom(IMG_DIR,GT_DIR,UCM_DIR,EV_DIR);
    eval_Fop(IMG_DIR, GT_DIR, UCM_DIR,EV_DIR);
    plot_eval(EV_DIR);
end
