% Evaluation script for crisp boundaries
%
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2015 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

%% paths (modify these to point where you want)
%% linux
githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir));cd(githubdir);
seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));% linux
DATA_DIR ='/home/lun5/HEproject/'; % linux
% %GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated_Oct', 'best_images_july30');
% %RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_PJoint_scale_offset_all3_newsetup');
% %RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_PMI_fullres_all3');
% %RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','oppcol_3channel_scale_offset_mult15');
IMG_DIR = fullfile(DATA_DIR,'data','normalization_512');
GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_fine_coarse');
RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','PMI_hue_lowres_accurate');
%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/Test';
IMG_DIR = '/home/lun5/HEproject/data/normalization_512/';
GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_fine_coarse');
RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','eval_PMI_hue_fullscale');

%% window
%githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; 
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
% IMG_DIR  = '/Users/lun5/Research/data/normalization_512';

% IMG_DIR = fullfile(DATA_DIR,'Tiles_512');
% IMG_DIR = fullfile(DATA_DIR,'Tiles_512','Test');
%GT_DIR = '/Users/lun5/Research/data/groundTruth/groundTruth_512_512';
% GT_DIR = '/Users/lun5/Research/data/groundTruth/groundTruth_512_512_reannotated_Oct/best_images_july30';
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','JSEG');
% bsrdir = '/Users/lun5/Research/packages/BSR/grouping';addpath(genpath(bsrdir));
%RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','PMI_hue_lowres_accurate');
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','PJoint_hue_fullres');
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','PJoint_hue_lowres');
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','PMI_hue_lowres_mult2');

%% set environment for affinity calculation
opts_affinity = setEnvironment_affinity;
opts_affinity.features.which_features = {'hue opp'};
%opts_affinity.features.which_features = {'hue opp', 'brightness opp', 'saturation opp'};
opts_affinity.joint_exponent = 2;%rho_list(i);
opts_affinity.sig = 3; %sigma_list(i); For 512x512 is 0.25, 2048 is 4
opts_affinity.display_progress = false;
opts_affinity.affinity.plot = false;
opts_affinity.features.plot = false;
opts_affinity.scale_offset = 1;
opts_affinity.affinityFunction = 'PMI';

opts_clustering = setEnvironment_clustering;
opts_clustering.display_progress = false;
opts_clustering.calculate_segments = false;
opts_clustering.plot_results = false;
opts_clustering.spectral_clustering.approximate = false;
%opts_affinity.scale_offset = 1;
%opts_affinity.affinityFunction = 'PMI';
evalAll(IMG_DIR,GT_DIR,RESULTS_DIR, opts_affinity, opts_clustering);
