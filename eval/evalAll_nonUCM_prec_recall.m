% evaluation script for regions based on Burak script

function [] = evalAll_nonUCM_prec_recall(IMG_DIR,GT_DIR,RESULTS_DIR,ev_mode)
    
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
    %img_list = dirrec(SEG_DIR,'.mat');
    %fprintf('Number of files is %d\n',length(img_list));
    EV_DIR = fullfile(RESULTS_DIR,['ev_txt_' ev_mode '_burak']);
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_all232');
    %EV_DIR = fullfile(RESULTS_DIR,'ev_txt_reannotated_Oct14');
    if (~exist(EV_DIR,'dir'))
        mkdir(EV_DIR);
    end
    [~,im_name,~] = fileparts(img_list{1});    
    tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    segs = tmp.data; 
    nSegments = length(segs); % segments 2:2:200
precision_mat = zeros(length(img_list), nSegments);
recall_mat = zeros(length(img_list), nSegments);
for i =1 :length(img_list)
    [~,im_name,~] = fileparts(img_list{i});    
    tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    segs = tmp.data; 
    nSegments = length(segs); % segments 2:2:200
    
    % load ground truth
    tmp = load(fullfile(GT_DIR_mode,[im_name '.mat']));
    groundTruth = tmp.groundTruth;
    for j = 1:nSegments
        [precision,recall] = evalPrecisionRecall(groundTruth,segs{j});
        precision_mat(i,j) = precision;
        recall_mat(i,j) = recall;
    end 
end

F_score_mat = 2*(precision_mat.*recall_mat)./((precision_mat+recall_mat) +...
    (sum(precision_mat(:)+recall_mat(:))==0));
F_score_avg = mean(F_score_mat,1);
[~,param_setting_ods] = max(F_score_avg);

recall_ods = recall_mat(:, param_setting_ods);
precision_ods =  precision_mat(:, param_setting_ods);
F_score_ods = F_score_mat(:,param_setting_ods);
%max
[F_score_ois, param_setting_ois] = max(F_score_mat,[],2);
recall_ois = recall_mat(sub2ind(size(recall_mat),1:length(img_list),param_setting_ois'));
precision_ois = precision_mat(sub2ind(size(recall_mat),1:length(img_list),param_setting_ois'));

f = fopen(fullfile(EV_DIR,'eval_summary.txt'),'w');
for i = 1 :length(img_list)
    [~,im_name,~] = fileparts(img_list{i});  
    fprintf(f,'%s %d %.2f %.2f %.2f %d %.2f %.2f %.2f\n',im_name,...
        param_setting_ods, precision_ods(i), recall_ods(i), F_score_ods(i),...
        param_setting_ois(i), precision_ois(i),recall_ois(i), F_score_ois(i));
end
fclose(f);

    %fprintf('nSegments is %d\n',nSegments);
    %nSegments = 1;
    %% eval using BSR metrics
