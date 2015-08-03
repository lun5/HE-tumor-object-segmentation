% generate output at OIS for best images in the test set
% Luong Nguyen 07/23/2015
% this is to test why the performanc isn't very good
DATA_DIR ='/home/lun5/HEproject/'; % linux
%IMG_DIR = fullfile(DATA_DIR,'data','images','test');
GT_DIR = fullfile(DATA_DIR,'data','groundTruth_512_512');
%IMG_DIR = '/home/lun5/HEproject/TilesForLabeling_tiff_renamed';
IMG_DIR = '/home/lun5/HEproject/data/images/test';
RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','combinedFeatures_32images');
gt_display = fullfile(DATA_DIR,'data','gt_1024_diplay');
outDir = fullfile(RESULTS_DIR,'outputTop');
if ~exist(outDir,'dir')
    mkdir(outDir)
end

eval_img = dlmread(fullfile(RESULTS_DIR,'ev_txt','eval_bdry_img.txt'));
%eval_img = dlmread(fullfile(RESULTS_DIR,'ev_txt','eval_Fop_img.txt'));
ind_top = 1:32;%find(eval_img(:,5) > 0.5);
best_thres = eval_img(ind_top,2);
IMG_EXT = '.tif';
img_list = dirrec(IMG_DIR,IMG_EXT);
img_list = img_list(:,ind_top);

parfor i = 1:numel(img_list)
    [~,im_name,~] = fileparts(img_list{i});
    weight_str = 'weights_3_1_0';
    outDir_wt = fullfile(outDir,[weight_str,'_bdry']);
    if ~ exist(outDir_wt,'dir')
        mkdir(outDir_wt)
    end
    fprintf('\n\nCalculate best segmentation for image %s...',im_name); T = tic;
    if ~exist(fullfile(outDir_wt,[im_name '_bestSeg.tif']),'file')
        tmp = load(fullfile(RESULTS_DIR,'ucm2',weight_str,[im_name '.mat']));
        ucm2 = tmp.data;
        I = imread(img_list{i});I = double(I(1:4:end,1:4:end,:));
        thr = best_thres(i);
        segmented_image = ucm2colorsegs(ucm2,I,thr);
        imwrite(uint8(segmented_image),fullfile(outDir_wt,[im_name '_bestSeg.tif']),'Resolution',300);
        J = label2rgb(bwlabel(ucm2(2:end,2:end) <= thr),'jet');
        imwrite(J,fullfile(outDir_wt,[im_name,'_bestSeg_colr.tif']),'Resolution',300);
    end
    
    filenames = {img_list{i},...
        fullfile(outDir_wt,[im_name '_bestSeg.tif']),...
        fullfile(outDir_wt,[im_name '_bestSeg_colr.tif']),...
        fullfile(gt_display,[lower(im_name) '.tif'])};
    a = cell(4,1);
    for j = 1:length(filenames)
        a{j} = imread(filenames{j});
    end
    a{1} = a{1}(1:2:end,1:2:end,:);
    montage = cat(2,a{:});
    imwrite(montage,fullfile(outDir_wt,[im_name '_montage.tif']),'Resolution',300);
    t = toc(T); fprintf('done: %1.2f sec\n', t);
end