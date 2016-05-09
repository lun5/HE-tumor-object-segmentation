% evaluation script for regions based on Burak script

function [] = evalAll_nonUCM_overlap(GT_DIR,RESULTS_DIR,ev_mode)
    
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
    SEG_DIR = fullfile(RESULTS_DIR,'segmented_images');
%     IMG_EXT = '.mat';
%     img_list = dirrec(SEG_DIR,'.mat');
    fprintf('Number of files is %d\n',length(img_list));
    EV_DIR = fullfile(RESULTS_DIR,['ev_txt_' ev_mode '_overlap_April4']);
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_all232');
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_reannotated_Oct14');
    if (~exist(EV_DIR,'dir'))
        mkdir(EV_DIR);
    end
    [~,im_name,~] = fileparts(img_list{1});    
    tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    %segs = tmp.data; 
    %nSegments = length(segs); % segments 2:2:200
    segs = tmp.data{1}; % Burak's results
    nSegments = length(segs); % segments 2:2:200
    ff_mat = zeros(length(img_list), nSegments);
    bb_mat = zeros(length(img_list), nSegments);
    fr_mat = zeros(length(img_list), nSegments);
    alpha = 0.75;
for i =1:length(img_list)
    [~,im_name,~] = fileparts(img_list{i}); 
    tic;
    tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    %segs = tmp.data; 
    %nSegments = length(segs); % segments 2:2:200
    segs = tmp.data{1}; % Burak's results
    nSegments = length(segs); % segments 2:2:200
    % load ground truth
    tmp = load(fullfile(GT_DIR_mode,[im_name '.mat']));
    groundTruth = tmp.groundTruth;
    parfor j = 1:nSegments
        seg = segs{j,1}; seg = uint8(seg(1:4:end,1:4:end))+1;%Burak
        %seg = segs{j}; 
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
% recall_pen_ods = recall_pen_mat(:, param_setting_pen_ods);
% precision_pen_ods =  precision_pen_mat(:, param_setting_pen_ods);
% F_score_pen_ods = F_score_pen_mat(:,param_setting_pen_ods);
% %max
% [F_score_pen_ois, param_setting_pen_ois] = max(F_score_pen_mat,[],2);
% recall_pen_ois = recall_pen_mat(sub2ind(size(recall_pen_mat),1:length(img_list),param_setting_pen_ois'));
% precision_pen_ois = precision_pen_mat(sub2ind(size(recall_pen_mat),1:length(img_list),param_setting_pen_ois'));
% 
% f_pen = fopen(fullfile(EV_DIR,'eval_pen_summary.txt'),'w');
% for i = 1 :length(img_list)
%     [~,im_name,~] = fileparts(img_list{i});  
%     fprintf(f_pen,'%s %d %.2f %.2f %.2f %d %.2f %.2f %.2f\n',im_name,...
%         param_setting_pen_ods, precision_pen_ods(i), recall_pen_ods(i), F_score_pen_ods(i),...
%         param_setting_pen_ois(i), precision_pen_ois(i),recall_pen_ois(i), F_score_pen_ois(i));
% end
% fclose(f_pen);

    %fprintf('nSegments is %d\n',nSegments);
    %nSegments = 1;
    %% eval using BSR metrics
