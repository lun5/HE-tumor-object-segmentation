%% script to evaluate the performance of gPb
% Luong Nguyen
% Adapted from BSR code, and Isola's crisp boundary
% 8/6/2015

%% paths (modify these to point where you want)
%{
%% linux
githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir));cd(githubdir);
seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));% linux
%bsrdir = '/home/lun5/github/BSR/grouping';addpath(genpath(bsrdir));
DATA_DIR ='/home/lun5/HEproject/'; % linux
%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/For_Om'; 
IMG_DIR = '/home/lun5/HEproject/data/normalization_512';
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated_Oct', 'best_images_july30');
GT_DIR = '/home/lun5/HEproject/groundTruth/coarse_fine_GT_512_512';
ev_mode = {'all_files','invasive','well_defined'};
%all_methods = {'MeanShift', 'EGB',fullfile('JSEG','scale1'), fullfile('JSEG','scale2'),fullfile('JSEG','scale3'),'Ncut','GraphRLM','cca_fused_white_purple','Gland_Seg'};
%all_methods = {'Isola_accurate_lowres_new'};
%all_methods = {'PMI_lowres_accurate'};
%all_methods ={'Isola_speedy_tight'};
%all_methods = {fullfile('JSEG_grayscale','scale1'), fullfile('JSEG_grayscale','scale2'),fullfile('JSEG_grayscale','scale3')};
%all_methods = {'updated_cca_voronoi'};
%all_methods = {'combined_1_1'};
% RESULTS_DIR = cell(length(all_methods),1);
% for i = 1:length(all_methods)
% 	%RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
% 	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
% end
% 
% for m = 1:length(ev_mode)
%   for i = 1 :length(all_methods)
%     if exist(fullfile(RESULTS_DIR{i},['ev_txt_' ev_mode{m}],'finalResults.txt'),'file')
%        continue;
%     end 
%     fprintf('Dir %s with ev_mode %s\n',RESULTS_DIR{i},ev_mode{m});
%     %evalAll_nonUCM(fullfile(IMG_DIR,ev_mode{m}),GT_DIR,RESULTS_DIR{i},ev_mode{m});
%     evalAll_bsr(fullfile(IMG_DIR,ev_mode{m}),GT_DIR,RESULTS_DIR{i},ev_mode{m});
%   end
% end
%all_methods = {'PMI_hue_20vec'};
%all_methods ={'Isola_speedy','eval_PMI_hue_offset','bsr'};
%all_methods = {'updated_cca_voronoi'};
%all_methods = {'SIC_1_standard','SIC_1_2_standard'};
%all_methods = {'MeanShift','ncut_multiscale_1_6',fullfile('JSEG','new_params','scale1'),fullfile('GraphRLM','new_params'),fullfile('EGB','seism_params')};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
	%RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end

for m = 3%1:length(ev_mode)
  for i = 1 :length(all_methods)
    if exist(fullfile(RESULTS_DIR{i},['ev_txt_' ev_mode{m}],'finalResults.txt'),'file')
       continue;
    end 
    fprintf('Dir %s with ev_mode %s\n',RESULTS_DIR{i},ev_mode{m});
    %evalAll_nonUCM(fullfile(IMG_DIR,ev_mode{m}),GT_DIR,RESULTS_DIR{i},ev_mode{m});
    evalAll_bsr(fullfile(IMG_DIR,ev_mode{m}),GT_DIR,RESULTS_DIR{i},ev_mode{m});
  end
end
%for i = 1:length(all_methods)
%   plot_eval(fullfile(RESULTS_DIR{i},'ev_txt'));
%end
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_non_expert','Maurice');
%evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR);

%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_non_expert','Om');
%evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR);
%}
%% window
githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; % window
seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism'; 
addpath(genpath(seismdir)); cd(githubdir);addpath(genpath(githubdir));
DATA_DIR = 'Z:\HEproject';
IMG_DIR = fullfile(DATA_DIR,'normalization_512');%'Z:\Tiles_512';
GT_DIR = fullfile(DATA_DIR,'data','GroundTruth','coarse_fine_GT_512_512');%Z:\HEproject\data\groundTruth_512_512';
%RESULTS_DIR = 'Z:\HEproject\evaluation_results\eval_non_expert\Om';
%all_methods = {'MeanShift', 'EGB',fullfile('JSEG','scale1'), ...
    %'ncut_multiscale_1_6','GraphRLM','GlandSeg','updated_cca_with_bg'};
%all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
%    'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg','SuperPixel',...
%    fullfile('eval_non_expert','Om'),fullfile('eval_non_expert','Maurice')};
%all_methods = {'GlandSeg','SuperPixel'};
all_methods = {'eval_PMI_hue_offset','Isola_speedy','bsr'};%,'PMI_lowres_accurate','SIC_1'};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	%RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end
ev_mode = {'all_files','well_defined','invasive'};

%evalAll_nonUCM_prec_recall(fullfile(IMG_DIR,ev_mode{3}),GT_DIR,RESULTS_DIR,ev_mode{3})
for m = 2:3
    for i =1:length(all_methods)
        fprintf('directory %s at mode %s ...',RESULTS_DIR{i},ev_mode{m});
        T = tic;
%         if exist(fullfile(RESULTS_DIR{i},['ev_txt_' ev_mode{m} '_burak'],'eval_summary.txt'),'file')
%             continue;
%         end
        %evalAll_nonUCM_prec_recall(fullfile(IMG_DIR,ev_mode{m}),GT_DIR,RESULTS_DIR{i},ev_mode{m});
        %evalAll_UCM_prec_recall(fullfile(IMG_DIR,ev_mode{m}),GT_DIR,RESULTS_DIR{i},ev_mode{m});
        %evalAll_nonUCM_overlap(GT_DIR,RESULTS_DIR{i},ev_mode{m});
        evalAll_UCM_overlap(GT_DIR,RESULTS_DIR{i},ev_mode{m});
        fprintf('Done in %.2f\n',toc(T));
    end
end
%% mac
% githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation'; % mac
% seismdir = '/Users/lun5/Research/github/seism'; %mac
% addpath(genpath(seismdir)); cd(githubdir)
% DATA_DIR = '/Users/lun5/Research/data';
% %IMG_DIR = fullfile(DATA_DIR,'Tiles_512');
% IMG_DIR = fullfile(DATA_DIR,'Tiles_512','Test');
% %GT_DIR = '/Users/lun5/Research/data/groundTruth/groundTruth_512_512';
% GT_DIR = '/Users/lun5/Research/data/groundTruth/groundTruth_512_512_reannotated/best_images_july30';
% %RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','JSEG');
% %bsrdir = '/Users/lun5/Research/packages/BSR/grouping';addpath(genpath(bsrdir));
% %evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR);
% evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR);
