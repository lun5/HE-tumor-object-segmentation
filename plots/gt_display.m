%% get ground truth display
gt_dir = '/Users/lun5/Research/data/groundTruth_1024_1024';
FILE_EXT = '.mat';
img_list = dirrec(gt_dir,FILE_EXT);
out_dir = '/Users/lun5/Research/data/groundTruth_1024_1024/diplay';
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end
    
parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        tmp = load(fullfile(gt_dir,[im_name '.mat']));
        seg = tmp.groundTruth{1}.Segmentation;
        imwrite(label2rgb(seg,'jet'),fullfile(out_dir,[im_name '.tif']),'Resolution',300);
end