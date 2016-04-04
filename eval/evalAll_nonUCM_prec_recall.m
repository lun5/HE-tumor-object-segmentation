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
%     IMG_EXT = '.mat';
%     img_list = dirrec(SEG_DIR,'.mat');
    fprintf('Number of files is %d\n',length(img_list));
    EV_DIR = fullfile(RESULTS_DIR,['ev_txt_' ev_mode '_burak_March30']);
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
    precision_pen_mat = zeros(length(img_list), nSegments);
    recall_pen_mat = zeros(length(img_list), nSegments);
    ri_mat = zeros(length(img_list), nSegments);
    gce_mat = zeros(length(img_list), nSegments);
    vi_mat =zeros(length(img_list), nSegments);
    pen_mat = zeros(length(img_list),nSegments);
for i =1:length(img_list)
    [~,im_name,~] = fileparts(img_list{i});    
    tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    segs = tmp.data; 
    nSegments = length(segs); % segments 2:2:200
    
    % load ground truth
    tmp = load(fullfile(GT_DIR_mode,[im_name '.mat']));
    groundTruth = tmp.groundTruth;
    for j = 1:nSegments
        %seg = segs{j,1}; seg = uint8(seg(1:4:end,1:4:end));
        seg = segs{j};
        [precision,recall, penalty, ri, gce, vi] = ...
             evalPrecisionRecall_new(groundTruth,seg);
        precision_mat(i,j) = precision;
        recall_mat(i,j) = recall;
        precision_pen_mat(i,j) = precision./max(penalty,1);
        recall_pen_mat(i,j) = recall./max(penalty,1);
        ri_mat(i,j) = ri;
        gce_mat(i,j) = gce;
        vi_mat(i,j) = vi;
        pen_mat(i,j) = penalty;
    end 
end

F_score_mat = 2*(precision_mat.*recall_mat)./((precision_mat+recall_mat) +...
    ((precision_mat+recall_mat)==0));
% F_score_avg = mean(F_score_mat,1);
% [~,param_setting_ods] = max(F_score_avg);
save(fullfile(EV_DIR,'ri_mat.mat'),'ri_mat');
save(fullfile(EV_DIR,'precision_mat.mat'),'precision_mat');
save(fullfile(EV_DIR,'recall_mat.mat'),'recall_mat');
save(fullfile(EV_DIR,'gce_mat.mat'),'gce_mat');
save(fullfile(EV_DIR,'vi_mat.mat'),'vi_mat');
save(fullfile(EV_DIR,'F_score_mat.mat'),'F_score_mat');
save(fullfile(EV_DIR,'pen_mat.mat'),'pen_mat');
save(fullfile(EV_DIR,'precision_pen_mat.mat'),'precision_pen_mat');
save(fullfile(EV_DIR,'recall_pen_mat.mat'),'recall_pen_mat');

ri_average = mean(ri_mat,1);
[~,param_setting_ods] = max(ri_average);
recall_ods = recall_mat(:, param_setting_ods);
precision_ods =  precision_mat(:, param_setting_ods);
F_score_ods = F_score_mat(:,param_setting_ods);

ri_ods = ri_mat(:,param_setting_ods);
gce_ods = gce_mat(:,param_setting_ods);
vi_ods = vi_mat(:,param_setting_ods);

%max
[F_score_ois, param_setting_ois] = max(F_score_mat,[],2);
recall_ois = recall_mat(sub2ind(size(recall_mat),1:length(img_list),param_setting_ois'));
precision_ois = precision_mat(sub2ind(size(precision_mat),1:length(img_list),param_setting_ois'));
ri_ois = ri_mat(sub2ind(size(precision_mat),1:length(img_list),param_setting_ois'));
gce_ois = gce_mat(sub2ind(size(precision_mat),1:length(img_list),param_setting_ois'));
vi_ois = vi_mat(sub2ind(size(precision_mat),1:length(img_list),param_setting_ois'));

f = fopen(fullfile(EV_DIR,'eval_summary_new.txt'),'w');
for i = 1 :length(img_list)
    [~,im_name,~] = fileparts(img_list{i});  
    fprintf(f,'%s %d %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f\n',im_name,...
        param_setting_ods, precision_ods(i), recall_ods(i), F_score_ods(i),...
        ri_ods(i), gce_ods(i), vi_ods(i),...
        param_setting_ois(i), precision_ois(i),recall_ois(i), F_score_ois(i),...
        ri_ois(i), gce_ois(i), vi_ois(i));
end
fclose(f);

F_score_pen_mat = 2*(precision_pen_mat.*recall_pen_mat)./((precision_pen_mat+recall_pen_mat) +...
    ((precision_pen_mat+recall_pen_mat)==0));
F_score_pen_avg = mean(F_score_pen_mat,1);
save(fullfile(EV_DIR,'F_score_pen_mat.mat'),'F_score_pen_mat');
[~,param_setting_pen_ods] = max(F_score_pen_avg);

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
