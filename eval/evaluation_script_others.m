%% script to evaluate the performance of gPb
% Luong Nguyen
% Adapted from BSR code, and Isola's crisp boundary
% 8/6/2015

%% paths (modify these to point where you want)
%% linux
githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir));cd(githubdir);
seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));% linux
bsrdir = '/home/lun5/github/BSR/grouping';addpath(genpath(bsrdir));
DATA_DIR ='/home/lun5/HEproject/'; % linux
%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/For_Om';
IMG_DIR = '/home/lun5/HEproject/data/normalization_512/Test';
GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated_Oct', 'best_images_july30');
%RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','Isola_speedy');
%RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','Isola_color_accurate_low_res');
%evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR);
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512');

all_methods = {'MeanShift', 'EGB',fullfile('JSEG','scale1'), fullfile('JSEG','scale2'),fullfile('JSEG','scale3'),'Ncut','connected_comp','GraphRLM'};
all_methods = {'cca_fused_white_purple'};
%all_methods = {'Gland_Seg'};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
end

%RESULTS_DIR{2} = fullfile(DATA_DIR,'evaluation_results','QuadTree');
%RESULTS_DIR{7} = fullfile(DATA_DIR,'evaluation_results','GraphRLM','new_params');
%RESULTS_DIR{8} = fullfile(DATA_DIR,'evaluation_results','bsr');
%RESULTS_DIR{9} = fullfile(DATA_DIR,'evaluation_results','Isola_multiscale');
%RESULTS_DIR{10} = fullfile(DATA_DIR,'evaluation_results','Isola_lowres_accurate');
%RESULTS_DIR{11} = fullfile(DATA_DIR,'evaluation_results','Isola_speedy');
%for i = 8:11
%    evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR{i});
%end
%
%for i = 1 :length(all_methods)
%    fprintf('Dir %s\n',RESULTS_DIR{i});
%    evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR{i});
%end

for i = 1:length(all_methods)
   plot_eval(fullfile(RESULTS_DIR{i},'ev_txt'));
end
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_non_expert','Maurice');
%evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR);

%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_non_expert','Om');
%evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR);

%% window
%githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; % window
%seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism'; 
%addpath(genpath(seismdir)); cd(githubdir)
% DATA_DIR = 'Z:\';
% IMG_DIR = 'Z:\Tiles_512';
% GT_DIR = 'Z:\HEproject\data\groundTruth_512_512';
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','EGB');

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
