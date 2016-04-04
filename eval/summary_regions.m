DATA_DIR = 'Z:\HEproject';
IMG_DIR = fullfile(DATA_DIR,'normalization_512');%'Z:\Tiles_512';
GT_DIR = fullfile(DATA_DIR,'data','GroundTruth','coarse_fine_GT_512_512');%Z:\HEproject\data\groundTruth_512_512';

all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg','SuperPixel',...
'eval_PMI_hue_offset','Isola_speedy','bsr'};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	%RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end
ev_mode = {'all_files','well_defined','invasive'};
% for i = 1:length(all_methods)
%     fprintf('directory %s ...\n',RESULTS_DIR{i});
%     T = tic;
%     load(fullfile(RESULTS_DIR{i},'ev_txt_well_defined_burak_April2','ri_mat.mat'))
%     well_defined = value_mat;
%     ri_average = mean(value_mat,1);
%     [~,param_setting_ods] = max(ri_average);
%     value_ods = value_mat(:,param_setting_ods);
%     fprintf('Well defined PRI score is %.4f\n',mean(value_ods));
%     load(fullfile(RESULTS_DIR{i},'ev_txt_invasive_burak_April2','ri_mat.mat'));
%     invasive_ri = value_mat;
%     ri_average = mean(value_mat,1);
%     [~,param_setting_ods] = max(ri_average);
%     value_ods = value_mat(:,param_setting_ods);
%     fprintf('invasive PRI score is %.4f\n',mean(value_ods));
%     value_mat = [well_defined; invasive_ri];
%     ri_average = mean(value_mat,1);
%     [~,param_setting_ods] = max(ri_average);
%     value_ods = value_mat(:,param_setting_ods);
%     fprintf('combined PRI score is %.4f\n',mean(value_ods));   
% end

for i = 1:7%length(all_methods)
    fprintf('directory %s ...\n',RESULTS_DIR{i});
    T = tic;
    %load(fullfile(RESULTS_DIR{i},'ev_txt_well_defined_burak_March30','precision_pen_mat.mat'))
    %load(fullfile(RESULTS_DIR{i},'ev_txt_well_defined_burak_March30','recall_pen_mat.mat'))
    %F_score_pen_mat = 2*(precision_pen_mat.*recall_pen_mat)./((precision_pen_mat+recall_pen_mat) +...
    %    ((precision_pen_mat+recall_pen_mat)==0));
    %load(fullfile(RESULTS_DIR{i},'ev_txt_well_defined_burak_March30','F_score_pen_mat.mat'));
    load(fullfile(RESULTS_DIR{i},'ev_txt_well_defined_overlap_April4','fr_mat.mat'));
    value_mat = fr_mat;%F_score_pen_mat;
    %load(fullfile(RESULTS_DIR{i},'ev_txt_well_defined_burak_April3','vi_mat.mat'));
    %value_mat = vi_mat;
    well_defined = value_mat;
    val_average = mean(value_mat,1);
    [~,param_setting_ods] = max(val_average);
    value_ods = value_mat(:,param_setting_ods);
    fprintf('Well defined score is %.4f\n',mean(value_ods));
    %load(fullfile(RESULTS_DIR{i},'ev_txt_invasive_burak_March30','precision_pen_mat.mat'));
    %load(fullfile(RESULTS_DIR{i},'ev_txt_invasive_burak_March30','recall_pen_mat.mat'));
    %F_score_pen_mat = 2*(precision_pen_mat.*recall_pen_mat)./((precision_pen_mat+recall_pen_mat) +...
    %    ((precision_pen_mat+recall_pen_mat)==0));
    %load(fullfile(RESULTS_DIR{i},'ev_txt_invasive_burak_March30','F_score_pen_mat.mat'));
    %value_mat = F_score_pen_mat;
    load(fullfile(RESULTS_DIR{i},'ev_txt_invasive_overlap_April4','fr_mat.mat'));
    value_mat = fr_mat;
    %load(fullfile(RESULTS_DIR{i},'ev_txt_invasive_burak_April3','vi_mat.mat'));
    %value_mat = vi_mat;
    invasive = value_mat;
    val_average = mean(value_mat,1);
    [~,param_setting_ods] = max(val_average);
    value_ods = value_mat(:,param_setting_ods);
    fprintf('invasive score is %.4f\n',mean(value_ods));
    value_mat = [well_defined; invasive];
    val_average = mean(value_mat,1);
    [~,param_setting_ods] = max(val_average);
    value_ods = value_mat(:,param_setting_ods);
    fprintf('combined score is %.4f\n',mean(value_ods));   
end