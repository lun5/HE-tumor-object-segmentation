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
    %img_list = dirrec(GT_DIR,'.mat');
    %% compute boundaries for images
    if (~exist(RESULTS_DIR,'dir'))
        mkdir(RESULTS_DIR);
    end
    
    SEG_DIR = fullfile(RESULTS_DIR,'segmented_images');
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_reannotated_Oct14');
    EV_DIR = fullfile(RESULTS_DIR,'ev_txt');
    if (~exist(EV_DIR,'dir'))
        mkdir(EV_DIR);
    end
    
    %nSegments = 100;% ncuts_color
    %nSegments = 200;% EGB
    %nSegments = 8;% JSEG multiscale
    %nSegments = 12;% JSEG one scale, GraphRLM

    [~,im_name,~] = fileparts(img_list{1});    
    tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    segs = tmp.data; 
    nSegments = length(segs); % segments 2:2:200

%     parfor i = 1:length(img_list)
%         [~,im_name,~] = fileparts(img_list{i});
%         fprintf('\n\nCalculate Ncuts Segmentation %s...',im_name);T = tic;
%         outFile = fullfile(RESULTS_DIR,'segmented_images',[im_name '.mat']);
%         if exist(outFile,'file')
%             continue;
%         end
%         I = double(imread(img_list{i}));
%         %I = mean(I,3); [nr,nc] = size(I);
%         seg = cell(nSegments/2,1);
%         [p,q,r] = size(I);
%         [layers,C]=compute_layers_C_multiscale(p,q);
%         dataW = computeParametersW(I);
%         %t =cputime;
%         W=computeMultiscaleW(I,layers,dataW,[]);
%         %fprintf('Calculate %d eigenvectors in %.2f seconds\n',nSegments,cputime-t);
%         if ~isempty(C)
%             [X,~,~] = computeNcutConstraint_projection(W,C,nSegments);
%         else
%             [X,~,~] = computeKFirstEigenvectors(W,nSegments);
%         end
%         indPixels = (1:p*q)';
%         X = reshape(X(indPixels,:),p,q,nSegments);
% 
%         %mult = 4; I = double(I(1:mult:end,1:mult:end,:));    
%         %[W,~] = ICgraph(I);
%         %[NcutEigenvectors,~] = ncut(W,nSegments + 1 + added_nsegs);
%         for j= 1:(nSegments/2)
%             %[NcutDiscrete,~,~] = ncutW(W,j);
%             %[NcutDiscrete,~] = discretisation(NcutEigenvectors(:,1:j + added_nsegs));
%             %SegLabel = zeros(nr,nc);
%             %for k=1:size(NcutDiscrete,2),
%             %    SegLabel = SegLabel + k*reshape(NcutDiscrete(:,k),nr,nc);
%             %end   
%             [SegLabel,~] =discretisation(X(:,:,1:2*j));                 
%             seg{j} = SegLabel;            
%         end
%         t = toc(T); fprintf('done: %1.2f sec\n', t);
%         parsave(outFile,seg);
%     end
%     
    %% eval using BSR metrics
    allBench_custom(IMG_DIR,GT_DIR,SEG_DIR,EV_DIR,nSegments);
    plot_eval(EV_DIR);


