%% get ground truth display
%gt_dir = '/Users/lun5/Research/data/groundTruth_1024_1024';
gt_dir = '/Users/lun5/Research/HE_Segmentation/groundTruth/coarse_fine_GT_512_512/all_files';
IM_GT_DIR = '/Users/lun5/Research/HE_Segmentation/groundTruth/coarse_fine_GT_512_512/bw_images';

FILE_EXT = '.mat';
img_list = dirrec(gt_dir,FILE_EXT);
%out_dir = '/Users/lun5/Research/data/groundTruth_1024_1024/diplay';
out_dir = '/Users/lun5/Research/HE_Segmentation/groundTruth/coarse_fine_GT_512_512/bw_images';
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end
se = strel('disk',2,4);

for i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        tmp = load(fullfile(gt_dir,[im_name '.mat']));
        seg = tmp.groundTruth{1}.Segmentation;
        bdry_comp = seg2bdry(seg,'imageSize');
        bdry_comp = repmat(uint8(imdilate(bdry_comp,se).*255),1,1,3);
        imwrite(bdry_comp,fullfile(out_dir,[im_name '.tif']));
        %imwrite(label2rgb(seg,'jet'),fullfile(out_dir,[im_name '.tif']),'Resolution',300);
end