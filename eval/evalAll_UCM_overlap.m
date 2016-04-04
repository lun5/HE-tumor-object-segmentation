% evaluation script for regions based on Burak script

function [] = evalAll_UCM_overlap(GT_DIR,RESULTS_DIR,ev_mode)
    
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
    SEG_DIR = fullfile(RESULTS_DIR,'ucm2');
%     IMG_EXT = '.mat';
%     img_list = dirrec(SEG_DIR,'.mat');
    fprintf('Number of files is %d\n',length(img_list));
    EV_DIR = fullfile(RESULTS_DIR,['ev_txt_' ev_mode '_overlap_April4']);
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_all232');
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_reannotated_Oct14');
    if (~exist(EV_DIR,'dir'))
        mkdir(EV_DIR);
    end
        
    nSegments = 50;
    fprintf('nSegments is %d\n',nSegments);
    thresh = linspace(1/(nSegments+1),1-1/(nSegments+1),nSegments)';
    save(fullfile(EV_DIR,'thresh'),'thresh');
    ff_mat = zeros(length(img_list), nSegments);
    bb_mat = zeros(length(img_list), nSegments);
    fr_mat = zeros(length(img_list), nSegments);
    
    alpha = 0.75;
for i =1:length(img_list)
    tic;
    [~,im_name,~] = fileparts(img_list{i});
    tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    ucm2 = tmp.data;
    ucm = double(ucm2);

    % load ground truth
    tmp = load(fullfile(GT_DIR_mode,[im_name '.mat']));
    groundTruth = tmp.groundTruth;
    parfor j = 1:nSegments
        labels2 = bwlabel(ucm <= thresh(j));
        seg = labels2(2:2:end, 2:2:end);
        [ff_score, bb_score] = evalRegions(groundTruth,seg);
        ff_mat(i,j) = ff_score;
        bb_mat(i,j) = bb_score;
        if bb_score == -1;
            fr_mat(i,j) = ff_score;
        else
            fr_mat(i,j) = alpha*ff_score + (1-alpha)*bb_score;
        end
    end 
    toc
end

save(fullfile(EV_DIR,'ff_mat.mat'),'ff_mat');
save(fullfile(EV_DIR,'bb_mat.mat'),'bb_mat');
save(fullfile(EV_DIR,'fr_mat.mat'),'fr_mat');

