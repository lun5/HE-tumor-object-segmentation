DATA_DIR = 'Z:\HEproject';
IMG_DIR = fullfile(DATA_DIR,'normalization_512');%'Z:\Tiles_512';
GT_DIR = fullfile(DATA_DIR,'data','GroundTruth','coarse_fine_GT_512_512');%Z:\HEproject\data\groundTruth_512_512';

all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'eval_PMI_hue_offset','Isola_speedy','bsr'};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	%RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end
ev_mode = {'all_files','well_defined','invasive'};
for i = 1:length(all_methods)
    fprintf('directory %s ...\n',RESULTS_DIR{i});
    T = tic;
    load(fullfile(RESULTS_DIR{i},'ev_txt_well_defined_burak','ri_mat.mat'))
    well_defined_ri = ri_mat;
    ri_average = mean(ri_mat,1);
    [~,param_setting_ods] = max(ri_average);
    ri_ods = ri_mat(:,param_setting_ods);
    fprintf('Well defined PRI score is %.4f\n',mean(ri_ods));
    load(fullfile(RESULTS_DIR{i},'ev_txt_invasive_burak','ri_mat.mat'));
    invasive_ri = ri_mat;
    ri_average = mean(ri_mat,1);
    [~,param_setting_ods] = max(ri_average);
    ri_ods = ri_mat(:,param_setting_ods);
    fprintf('invasive PRI score is %.4f\n',mean(ri_ods));
    ri_mat = [well_defined_ri; invasive_ri];
    ri_average = mean(ri_mat,1);
    [~,param_setting_ods] = max(ri_average);
    ri_ods = ri_mat(:,param_setting_ods);
    fprintf('combined PRI score is %.4f\n',mean(ri_ods));   
end