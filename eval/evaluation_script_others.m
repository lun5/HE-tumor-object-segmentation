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
%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/Test';
IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/';
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated_Oct', 'best_images_july30');
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512');
GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_fine_coarse');
RESULTS_DIR{1} = fullfile(DATA_DIR,'evaluation_results','MeanShift');
RESULTS_DIR{2} = fullfile(DATA_DIR,'evaluation_results','QuadTree');
RESULTS_DIR{3} = fullfile(DATA_DIR,'evaluation_results','EGB','seism_params');
RESULTS_DIR{4} = fullfile(DATA_DIR,'evaluation_results','JSEG','new_params','scale1');
RESULTS_DIR{5} = fullfile(DATA_DIR,'evaluation_results','JSEG','new_params','scale2');
RESULTS_DIR{6} = fullfile(DATA_DIR,'evaluation_results','ncut_multiscale_1_6');
RESULTS_DIR{7} = fullfile(DATA_DIR,'evaluation_results','GraphRLM','new_params');
RESULTS_DIR{8} = fullfile(DATA_DIR,'evaluation_results','bsr');
RESULTS_DIR{9} = fullfile(DATA_DIR,'evaluation_results','Isola_multiscale');
RESULTS_DIR{10} = fullfile(DATA_DIR,'evaluation_results','Isola_lowres_accurate');
RESULTS_DIR{11} = fullfile(DATA_DIR,'evaluation_results','Isola_speedy');

%evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR);
for i = 1:7
    evalAll_nonUCM(IMG_DIR,GT_DIR,RESULTS_DIR{i});
end

for i = 8:11
    evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR{i});
end


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

