%% script to evaluate the performance of gPb
% Luong Nguyen
% Adapted from BSR code, and Isola's crisp boundary
% 8/6/2015

%% paths (modify these to point where you want)
%githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; % window
%addpath(genpath(githubdir));
%seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));
seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism'; 
addpath(genpath(seismdir)); cd(githubdir)
%bsrdir = '/home/lun5/github/BSR/grouping';addpath(genpath(bsrdir));
DATA_DIR ='/home/lun5/HEproject/'; % linux
IMG_DIR = '/home/lun5/HEproject/data/Tiles_512';
%IMG_DIR = '/home/lun5/HEproject/data/images/test';
GT_DIR = fullfile(DATA_DIR,'data','groundTruth_512_512');
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','ncuts');
% DATA_DIR = 'Z:\';
% IMG_DIR = 'Z:\Tiles_512';
% GT_DIR = 'Z:\HEproject\data\groundTruth_512_512';
% RESULTS_DIR = 'Z:\HEproject\evaluation_results\ncuts_color';
RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','EGB');
%evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR);
evalAll_ncuts(IMG_DIR,GT_DIR,RESULTS_DIR);
