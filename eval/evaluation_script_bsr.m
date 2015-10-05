%% script to evaluate the performance of gPb
% Luong Nguyen
% Adapted from BSR code, and Isola's crisp boundary
% 8/6/2015

%% paths (modify these to point where you want)
%% linux
%githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
%addpath(genpath(githubdir));
%seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));% linux
%addpath(genpath(seismdir)); cd(githubdir)
%DATA_DIR ='/home/lun5/HEproject/'; % linux
%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512';
%GT_DIR = fullfile(DATA_DIR,'data','groundTruth_512_512');
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','ncuts');

%% window
%githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; % window
%seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism'; 
%addpath(genpath(seismdir)); cd(githubdir)
% DATA_DIR = 'Z:\';
% IMG_DIR = 'Z:\Tiles_512';
% GT_DIR = 'Z:\HEproject\data\groundTruth_512_512';
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','EGB');
copyfile(
%% mac
githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation'; % mac
seismdir = '/Users/lun5/Research/github/seism'; %mac
addpath(genpath(seismdir)); cd(githubdir)
DATA_DIR = '/Users/lun5/Research/data';
IMG_DIR = fullfile(DATA_DIR,'Tiles_512');
%bsrdir = '/home/lun5/github/BSR/grouping';addpath(genpath(bsrdir));
GT_DIR = '/Users/lun5/Research/data/groundTruth/groundTruth_512_512';
% RESULTS_DIR = 'Z:\HEproject\evaluation_results\ncuts_color';
RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','EGB');
%evalAll_bsr(IMG_DIR,GT_DIR,RESULTS_DIR);
evalAll_ncuts(IMG_DIR,GT_DIR,RESULTS_DIR);
