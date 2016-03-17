gt_dir = 'Z:\HEproject\data\GroundTruth\coarse_fine_GT_512_512_new\well_defined';
im_name = '9uixINHtjjiS';
result_dir = 'Z:\HEproject\object_proposals\updated_cca_with_bg\segmented_images';

tmp = load(fullfile(gt_dir,[im_name '.mat']));
groundTruth = tmp.groundTruth;

tmp = load(fullfile(result_dir,[im_name '.mat']));
seg = tmp.data;
[precision,recall] = evalPrecisionRecall(groundTruth,seg{1});