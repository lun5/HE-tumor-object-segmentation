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
IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/Test';
GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated_Oct', 'best_images_july30');
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_PJoint_scale_offset_all3_newsetup');
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_PMI_fullres_all3');
%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/';
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512');
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
% %RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','EGB');
% %RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','JSEG');
% %bsrdir = '/Users/lun5/Research/packages/BSR/grouping';addpath(genpath(bsrdir));
% %evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR);
% %evalAll_ncuts(IMG_DIR,GT_DIR,RESULTS_DIR);
RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_PMI_scale_offset_all3_newsetup');

%% set environment for affinity calculation
opts_affinity = setEnvironment_affinity;
%opts_affinity.features.which_features = {'hue opp'};
opts_affinity.features.which_features = {'hue opp', 'brightness opp', 'saturation opp'};
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','combinedFeatures_32images');
%% this still can be changed to get the best sigma, might try 0.25
opts_affinity.joint_exponent = 2;%rho_list(i);
opts_affinity.sig = 3; %sigma_list(i); For 512x512 is 0.25, 2048 is 4
opts_affinity.display_progress = false;
opts_affinity.affinity.plot = false;
opts_affinity.features.plot = false;
opts_affinity.scale_offset = 1;
%opts_affinity.affinityFunction = 'PJoint';
%evalAll(IMG_DIR,GT_DIR,RESULTS_DIR, opts_affinity);
evalAll_new(IMG_DIR,GT_DIR,RESULTS_DIR, opts_affinity);

%% Below are the benchmark numbers you should get for each type of parameter settings
%{
../Results/speedy
Boundary
ODS: F( 0.74, 0.72 ) = 0.73   [th = 0.07]
OIS: F( 0.76, 0.75 ) = 0.75
Area_PR = 0.78

Region
GT covering: ODS = 0.59 [th = 0.13]. OIS = 0.64. Best = 0.74
Rand Index: ODS = 0.82 [th = 0.10]. OIS = 0.85.
Var. Info.: ODS = 1.70 [th = 0.19]. OIS = 1.51.



../Results/accurate_low_res
Boundary
ODS: F( 0.73, 0.75 ) = 0.74   [th = 0.08]
OIS: F( 0.76, 0.77 ) = 0.77
Area_PR = 0.80

Region
GT covering: ODS = 0.61 [th = 0.13]. OIS = 0.66. Best = 0.76
Rand Index: ODS = 0.83 [th = 0.13]. OIS = 0.86.
Var. Info.: ODS = 1.58 [th = 0.20]. OIS = 1.42.



../Results/accurate_high_res
Boundary
ODS: F( 0.73, 0.72 ) = 0.73   [th = 0.06]
OIS: F( 0.77, 0.73 ) = 0.75
Area_PR = 0.78

Region
GT covering: ODS = 0.58 [th = 0.08]. OIS = 0.63. Best = 0.72
Rand Index: ODS = 0.81 [th = 0.04]. OIS = 0.85.
Var. Info.: ODS = 1.75 [th = 0.11]. OIS = 1.54.



../Results/accurate_multiscale
Boundary
ODS: F( 0.72, 0.73 ) = 0.73   [th = 0.04]
OIS: F( 0.77, 0.76 ) = 0.76
Area_PR = 0.77

Region
GT covering: ODS = 0.60 [th = 0.05]. OIS = 0.65. Best = 0.74
Rand Index: ODS = 0.82 [th = 0.04]. OIS = 0.85.
Var. Info.: ODS = 1.68 [th = 0.07]. OIS = 1.47.



../Results/MS_algorithm_from_paper
Boundary
ODS: F( 0.75, 0.72 ) = 0.73   [th = 0.06]
OIS: F( 0.77, 0.76 ) = 0.77
Area_PR = 0.78

Region
GT covering: ODS = 0.60 [th = 0.09]. OIS = 0.66. Best = 0.76
Rand Index: ODS = 0.82 [th = 0.06]. OIS = 0.86.
Var. Info.: ODS = 1.65 [th = 0.13]. OIS = 1.43.


%}
