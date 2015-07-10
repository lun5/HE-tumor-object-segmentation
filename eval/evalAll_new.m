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
        mkdir(fullfile(RESULTS_DIR,'E_oriented'));
        mkdir(fullfile(RESULTS_DIR,'ucm2'));
        mkdir(fullfile(RESULTS_DIR,'segmented_images'));
        mkdir(fullfile(RESULTS_DIR,'edgemap'));
        mkdir(fullfile(RESULTS_DIR,'montageImageSegment'));
    end
    %opts_affinity = setEnvironment_affinity;
    opts_clustering = setEnvironment_clustering;
    opts_clustering.display_progress = false;
    %opts_clustering.calculate_segments = false;
    opts_clustering.plot_results = false;  
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        if (~exist(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),'file'))
            fprintf('\n\nCalculate E oriented, E_ucm, segmented_image %s...',im_name); tic;
            I = imread(img_list{i});
            mult = 4; dz_im = I(1:mult:end,1:mult:end,:);
            I = double(dz_im);
            %[A,im_sizes] = getW(I,opts_affinity);[~, ~, E_oriented] = graphSegmentation(A,im_sizes,I,opts_clustering);
            [Ws,Ws_each_feature_set, im_sizes] = getW(I,opts_affinity);
            [segmented_image,E_ucm, E, E_oriented] = graphSegmentation(Ws,Ws_each_feature_set,im_sizes,I,opts)            
            %E_oriented = imresize(E_oriented,size(I(:,:,1)));
            %E = max(E_oriented,[],3); 
            imwrite(mat2gray(1-E),fullfile(RESULTS_DIR,'edgemap',[im_name,'_edgemap.tif']),'Resolution',300); 
            parsave(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),E_oriented);
            parsave(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']),E_ucm);
            imwrite(uint8(segmented_image),fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif']),'Resolution',300); 
            t = toc; fprintf('done: %1.2f sec\n', t);
        end      
    end
    
    %% need to add the normalizing step
%     %% normalize output scale
%     E_orienteds = cell(1,length(img_list));
%     parfor i=1:length(img_list)
%         [~,im_name,~] = fileparts(img_list{i});
%         tmp = load(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']));
%         E_orienteds{i} = tmp.data;
%     end
%     %
%     max_val = max(cellfun(@max_all,E_orienteds));
%     parfor i=1:length(img_list)
%         E_orienteds{i} = E_orienteds{i}/max_val;
%     end
%     %% run UCM on boundary maps
%     parfor i=1:length(img_list)
%         [~,im_name,~] = fileparts(img_list{i});
%         if (~exist(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']),'file'))
%             fprintf('\n\nCalculate UCM %s...',im_name); tic;
%             ucm2 = contours2ucm_crisp_boundaries(mat2gray(E_orienteds{i}),'doubleSize');
%             parsave(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']),ucm2);
%         end
%         if ~exist(fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif']),'file')
%             fprintf('\n\nCalculate segmented image %s...',im_name); tic;
%             tmp = load(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']));
%             ucm2 = tmp.data;
%             I = imread(img_list{i});I = double(I(1:4:end,1:4:end,:));
%             segmented_image = ucm2colorsegs(ucm2,I,0.2);
%             imwrite(uint8(segmented_image),fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif']),'Resolution',300); 
%             t = toc; fprintf('done: %1.2f sec\n', t);
%         end
%     end
    
    %% eval using BSR metrics
    allBench_custom(IMG_DIR,GT_DIR,fullfile(RESULTS_DIR,'ucm2'),RESULTS_DIR);
    plot_eval(RESULTS_DIR);
end
