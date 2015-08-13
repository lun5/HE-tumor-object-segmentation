% Evaluation script for crisp boundaries
%
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2015 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

%% setup
% first cd to the 'Code' directory, then run:
%compile; % this will check to make sure everything is compiled properly; if it is not, will try to compile it

%% choose parameter settings to evaluate
%type = 'speedy'; 
%type = 'accurate_low_res';
%type = 'accurate_high_res';
%type = 'accurate_multiscale';
%type = 'MS_algorithm_from_paper';

%% paths (modify these to point where you want)
%githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir));
%seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));
seismdir = '/Users/lun5/Research/github/seism'; addpath(genpath(seismdir));
cd(githubdir)
%DATA_DIR = fullfile(pwd,'data'); %'PATH/TO/BSDS';
%DATA_DIR ='/home/lun5/HEproject/'; % linux
DATA_DIR = '/Users/lun5/Research/data/';
%IMG_DIR = fullfile(DATA_DIR,'data','images');
%IMG_DIR = '/home/lun5/HEproject/TilesForLabeling_tiff_renamed';
IMG_DIR = fullfile(DATA_DIR,'TilesForLabeling_tiff_renamed','test');%'/home/lun5/HEproject/data/images/test';
GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated','best_images_july30');%fullfile(DATA_DIR,'data','groundTruth_512_512');
RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_reannotated');
%RESULTS_DIR = fullfile(pwd,'results','eval_col_val');
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_hue_512_512'); 
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_hue_2048_2048'); 
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_hue_512_512_sig3_exp2_newAffinity_Merge'); 
%%
%joint_exponent_list = [1 1.25 1.5 2 3];
%sigma_sample_dist_list = [0.25 0.5 1 2 3 5 10 15];

%[p, q] = meshgrid(joint_exponent_list, sigma_sample_dist_list);
%rho_list = p(:); sigma_list = q(:);
%numCombs = length(rho_list);
%parfor i = 2:numCombs
    %% set environment for affinity calculation
    opts_affinity = setEnvironment_affinity;
    %opts_affinity.features.which_features = {'luminance'};
    opts_affinity.which_features = {'hue opp', 'brightness opp', 'saturation opp'}; 
    %RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','combinedFeatures_32images');
    opts_affinity.joint_exponent = 2;%rho_list(i);
    opts_affinity.sig = 3; %sigma_list(i); For 512x512 is 0.25, 2048 is 4
    opts_affinity.display_progress = false;
    opts_affinity.affinity.plot = false;
    opts.features.plot = false;
    %evalAll(IMG_DIR,GT_DIR,RESULTS_DIR, opts_affinity);
    evalAll_new(IMG_DIR,GT_DIR,RESULTS_DIR, opts_affinity);
    %% clean up no need since rad and rad_inner does not change
    %delete(sprintf('%s/caches/ii_jj_caches/512_512.mat', pwd));
%end

% DATA_DIR = '/Users/lun5/Research/data/';
% %IMG_DIR = fullfile(DATA_DIR,'data','images');
% %IMG_DIR = '/home/lun5/HEproject/TilesForLabeling_tiff_renamed';
% IMG_DIR = fullfile(DATA_DIR,'TilesForLabeling_tiff_renamed');%'/home/lun5/HEproject/data/images/test';
% GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated','best_images_july30');%fullfile(DATA_DIR,'data','groundTruth_512_512');
% 
% iids = dir(fullfile(GT_DIR,'*.mat'));
% test_img_dir = fullfile(IMG_DIR,'test');
% if ~exist(test_img_dir,'dir'), mkdir(test_img_dir); end;
% parfor i = 1:length(iids)
%     fname = fullfile(IMG_DIR, strcat(iids(i).name(1:end-4), '.tif'));
%     copyfile(fname,fullfile(test_img_dir, strcat(iids(i).name(1:end-4), '.tif')));
% end
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
