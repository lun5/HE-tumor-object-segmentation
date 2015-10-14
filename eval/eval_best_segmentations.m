% generate output at OIS for best images in the test set
% Luong Nguyen 07/23/2015
% this is to test why the performanc isn't very good
% UPDATE: 10/12/15

%% mac
% githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation';
% addpath(genpath(githubdir)); cd(githubdir)
% seismdir = '/Users/lun5/Research/github/seism'; addpath(genpath(seismdir));
% DATA_DIR = '/Users/lun5/Research/data/';
% IMG_DIR = fullfile(DATA_DIR,'TilesForLabeling_tiff_renamed','test');%'/home/lun5/HEproject/data/images/test';
% GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated','best_images_july30');%fullfile(DATA_DIR,'data','groundTruth_512_512');
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_reannotated');
%% window
githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; % window
seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism'; 
addpath(genpath(seismdir)); cd(githubdir)
DATA_DIR = 'Z:\';
IMG_DIR = 'Z:\Tiles_512\Test';
GT_DIR = 'Z:\HEproject\data\groundTruth_512_512';
RESULTS_DIR = fullfile(DATA_DIR,'HEproject','evaluation_results','Isola_lowres_accurate');
gt_display = fullfile(DATA_DIR,'groundTruth','groundTruth_reannotated_display');
outDir = fullfile(RESULTS_DIR,'best_boundary');

if ~exist(outDir,'dir')
    mkdir(outDir)
end

eval_img = dlmread(fullfile(RESULTS_DIR,'ev_txt_reannotated','eval_bdry_img.txt'));
%eval_img = dlmread(fullfile(RESULTS_DIR,'ev_txt','eval_Fop_img.txt'));
best_thres = eval_img(:,2);
IMG_EXT = '.tif';
img_list = dirrec(IMG_DIR,IMG_EXT);
if ~ exist(outDir,'dir')
    mkdir(outDir)
end

parfor i = 1:numel(img_list)
    [~,im_name,~] = fileparts(img_list{i});   
    fprintf('\n\nCalculate best segmentation for image %s...',im_name); T = tic;
    outFile = fullfile(outDir,[im_name, '.jpg']);
    if ~exist(outFile,'file')
        I = imread(img_list{i});
        if exist(fullfile(RESULTS_DIR,'ucm2'),'dir')
            tmp = load(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']));
            ucm2 = tmp.data; ucm2 = ucm2(3:2:end,3:2:end);
            thr = best_thres;
            edge_map = (ucm2>=thr);
        else
            tmp = load(fullfile(RESULTS_DIR,'segmented_images',[im_name '.mat']));
            segs = tmp.data;
            thr = round(best_thres);
            edge_map = edge(segs{thr})
        end      
                edge_map = imdilate(edge_map, strel('disk',1));
        edge_map_im = I.*uint8(repmat(~edge_map,[1 1 3]));
        imwrite(edge_map_im,outFile);
        %imwrite(label2rgb(labels),fullfile(output_dir,'seg_im',[im_name, '_' num2str(q_thresh), '_seg.jpg']));
    end    
    t = toc(T); fprintf('done: %1.2f sec\n', t);
end