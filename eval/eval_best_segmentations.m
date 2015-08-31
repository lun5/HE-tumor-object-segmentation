% generate output at OIS for best images in the test set
% Luong Nguyen 07/23/2015
% this is to test why the performanc isn't very good
% DATA_DIR ='/home/lun5/HEproject/'; % linux
% %IMG_DIR = fullfile(DATA_DIR,'data','images','test');
% GT_DIR = fullfile(DATA_DIR,'data','groundTruth_512_512');
% %IMG_DIR = '/home/lun5/HEproject/TilesForLabeling_tiff_renamed';
% IMG_DIR = '/home/lun5/HEproject/data/images/test';
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','combinedFeatures_32images');

githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir));
seismdir = '/Users/lun5/Research/github/seism'; addpath(genpath(seismdir));
cd(githubdir)
DATA_DIR = '/Users/lun5/Research/data/';
IMG_DIR = fullfile(DATA_DIR,'TilesForLabeling_tiff_renamed','test');%'/home/lun5/HEproject/data/images/test';
GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated','best_images_july30');%fullfile(DATA_DIR,'data','groundTruth_512_512');
RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_reannotated');
gt_display_old = fullfile(DATA_DIR,'groundTruth','groundTruth_512_display');
gt_display_new = fullfile(DATA_DIR,'groundTruth','groundTruth_reannotated_display');
outDir = fullfile(RESULTS_DIR,'outputTop_015');
if ~exist(outDir,'dir')
    mkdir(outDir)
end

eval_img = dlmread(fullfile(RESULTS_DIR,'ev_txt','eval_bdry_img.txt'));
%eval_img = dlmread(fullfile(RESULTS_DIR,'ev_txt','eval_Fop_img.txt'));
ind_top = 1:41;%find(eval_img(:,5) > 0.5);
best_thres = eval_img(ind_top,2);
IMG_EXT = '.tif';
img_list = dirrec(IMG_DIR,IMG_EXT);
img_list = img_list(:,ind_top);
weight_str = 'weights_3_1_1';
outDir_wt = fullfile(outDir,[weight_str,'_bdry']);
if ~ exist(outDir_wt,'dir')
    mkdir(outDir_wt)
end
parfor i = 1:numel(img_list)
    [~,im_name,~] = fileparts(img_list{i});   
    fprintf('\n\nCalculate best segmentation for image %s...',im_name); T = tic;
    if ~exist(fullfile(outDir_wt,[im_name '_bestSeg.tif']),'file')
        tmp = load(fullfile(RESULTS_DIR,'ucm2',weight_str,[im_name '.mat']));
        ucm2 = tmp.data;ucm2 = ucm2(3:2:end,3:2:end);
        I = imread(img_list{i});I = double(I(1:4:end,1:4:end,:));
        thr = 0.15; %best_thres(i);
        segmented_image = ucm2colorsegs(ucm2,I,thr);
        imwrite(uint8(segmented_image),fullfile(outDir_wt,[im_name '_bestSeg.tif']));%,'Resolution',300);
        %J = label2rgb(bwlabel(ucm2(2:end,2:end) <= thr),'jet');
        %imwrite(J,fullfile(outDir_wt,[im_name,'_bestSeg_colr.tif']),'Resolution',300);
    end
    
    filenames = {img_list{i},...
        fullfile(outDir_wt,[im_name '_bestSeg.tif']),...
        fullfile(gt_display_old,[lower(im_name) '.bmp']),...
        fullfile(gt_display_new,[lower(im_name) '.bmp'])};
        %fullfile(outDir_wt,[im_name '_bestSeg_colr.tif']),...
        %fullfile(gt_display_old,[lower(im_name) '.tif'])};
    a = cell(4,1);
    for j = 1:length(filenames)
        a{j} = imread(filenames{j});
    end
    a{1} = a{1}(1:4:end,1:4:end,:);a{2} = a{2}(1:2:end, 1:2:end,:);
    montage = cat(2,a{:});
    imwrite(montage,fullfile(outDir_wt,[im_name '_montage.tif']));%,'Resolution',300);
    t = toc(T); fprintf('done: %1.2f sec\n', t);
end