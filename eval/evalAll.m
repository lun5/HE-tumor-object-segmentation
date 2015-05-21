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


function [] = evalAll(IMG_DIR,GT_DIR,RESULTS_DIR)
    
    %% read images
    IMG_EXT = '.tif';
    img_list = dirrec(IMG_DIR,IMG_EXT);

    %% compute boundaries for images
    if (~exist(RESULTS_DIR,'dir'))
        mkdir(RESULTS_DIR);
    end
    opts_affinity = setEnvironment_affinity;
    opts_clustering = setEnvironment_clustering;
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        if (~exist(fullfile(RESULTS_DIR,[im_name '_E_oriented.mat']),'file'))
            fprintf('\n\nCalculate E oriented %s...',im_name); tic;
            I = imread(img_list{i});dz_im = I(1:4:end,1:4:end,:);
            I = double(dz_im);[A,im_sizes] = getW(I,opts_affinity);
            [~, ~, E_oriented] = graphSegmentation(A,im_sizes,I,opts_clustering);
            E_oriented = imresize(E_oriented,size(I(:,:,1)));
            parsave(fullfile(RESULTS_DIR,[im_name '_E_oriented.mat']),E_oriented);
            t = toc; fprintf('done: %1.2f sec\n', t);
        end      
    end
    %% normalize output scale
    E_orienteds = cell(1,length(img_list));
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        tmp = load(fullfile(RESULTS_DIR,[im_name '_E_oriented.mat']));
        E_orienteds{i} = tmp.data;
    end
    %
    max_val = max(cellfun(@max_all,E_orienteds));
    parfor i=1:length(img_list)
        E_orienteds{i} = E_orienteds{i}/max_val;
    end
    %% run UCM on boundary maps
    parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        if (~exist(fullfile(RESULTS_DIR,[im_name '.mat']),'file'))
            fprintf('\n\nCalculate UCM %s...',imname); tic;
            ucm2 = contours2ucm_crisp_boundaries(E_orienteds{i},'doubleSize');
            parsave(fullfile(RESULTS_DIR,[im_name '.mat']),'ucm2');
            t = toc; fprintf('done: %1.2f sec\n', t);
        end
    end
    
    %% eval using BSR metrics
    allBench_custom(IMG_DIR,GT_DIR,RESULTS_DIR,RESULTS_DIR);
    plot_eval(RESULTS_DIR);
end