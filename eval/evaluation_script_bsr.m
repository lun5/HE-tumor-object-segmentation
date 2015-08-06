%% script to evaluate the performance of gPb
% Luong Nguyen
% Adapted from BSR code, and Isola's crisp boundary
% 8/6/2015

%% paths (modify these to point where you want)
githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir));
seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));
cd(githubdir)
bsrdir = '/home/lun5/github/BSR/grouping';
addpath(genpath(bsrdir));
DATA_DIR ='/home/lun5/HEproject/'; % linux
IMG_DIR = '/home/lun5/HEproject/data/Tiles_512';
%IMG_DIR = '/home/lun5/HEproject/data/images/test';
GT_DIR = fullfile(DATA_DIR,'data','groundTruth_512_512');
RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','bsr');
evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR);

