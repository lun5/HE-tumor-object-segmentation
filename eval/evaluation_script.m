% Evaluation script for crisp boundaries
%
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2015 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

%% paths (modify these to point where you want)
%% linux
%githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
%addpath(genpath(githubdir));cd(githubdir);
%seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));% linux
%DATA_DIR ='/home/lun5/HEproject/'; % linux
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_PJoint_scale_offset_all3_newsetup');
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_PMI_fullres_all3');
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_fine_coarse');
% %RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','oppcol_3channel_scale_offset_mult15');
%IMG_DIR = fullfile(DATA_DIR,'data','normalization_2048');
%IMG_DIR = fullfile(DATA_DIR,'data','normalization_512','all_files');
%RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','PMI_fullres_speedy');
%RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','colorStats_scan');
%GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','all_files'); 
%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512';
%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/Test';
%IMG_DIR = '/home/lun5/HEproject/data/normalization_512/';

%% window
% githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; 
% seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism'; 
% addpath(genpath(seismdir)); cd(githubdir)
% DATA_DIR = 'D:\Documents\HE_Segmentation';
% IMG_DIR = fullfile(DATA_DIR,'data','normlization_512');
% GT_DIR = 'Z:\HEproject\data\groundTruth_512_512';
% RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','colorStats_param_scan');
% 
% if ~exist(RESULTS_DIR,'dir')
%     mkdir(RESULTS_DIR);
% end

%% mac
githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation'; % mac
seismdir = '/Users/lun5/Research/github/seism'; %mac
addpath(genpath(seismdir)); cd(githubdir)
DATA_DIR = '/Users/lun5/Research/HE_Segmentation/';
IMG_DIR  = fullfile(DATA_DIR,'normalization_512');
GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','all_files');
RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','colorStats_param_scan');
if ~exist(RESULTS_DIR,'dir')
    mkdir(RESULTS_DIR);
end

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
%opts_affinity.joint_exponent = 2;%rho_list(i);
%opts_affinity.sig = 3; %sigma_list(i); For 512x512 is 0.25, 2048 is 4
opts_affinity.display_progress = false;
opts_affinity.affinity.plot = false;
opts_affinity.features.plot = false;
opts_affinity.scale_offset = 1;
opts_affinity.affinityFunction = 'PMI';

opts_clustering = setEnvironment_clustering;
opts_clustering.display_progress = false;
opts_clustering.calculate_segments = false;
opts_clustering.plot_results = false;
opts_clustering.spectral_clustering.approximate = true;
opts_clustering.spectral_clustering.nvec = 100;
%opts_affinity.scale_offset = 1;
%opts_affinity.affinityFunction = 'PMI';

joint_exponent_vec = [1.25, 2, 2.5, 3];
sig_vec = [0.25, 1, 3, 5, 7];

run_times = zeros(length(joint_exponent_vec)*length(sig_vec),1);
count = 0;
for i = 1:length(joint_exponent_vec)    
    for j = 1:length(sig_vec)
        opts_affinity.joint_exponent = joint_exponent_vec(i);
        opts_affinity.sig = sig_vec(j); 
        output_dir = fullfile(RESULTS_DIR, ['exp' '_' num2str(joint_exponent_vec(i))...
            '_sig_' num2str(sig_vec(j))]);
        tic;
        evalAll(IMG_DIR,GT_DIR,output_dir, opts_affinity, opts_clustering);
        count = count + 1;
        run_times(count) = toc;
        fprintf('Done with combo joint exp = %.2f, sig = %.2f in %2.f seconds\n',...
            joint_exponent_vec(i), sig_vec(j), run_times(count));
    end    
end
