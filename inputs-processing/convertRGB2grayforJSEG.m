% make grayscale images
im_dir = 'Z:\HEproject\data\normalization_512_jpg\';
im_list = dir(fullfile(im_dir,'*.jpg'));
im_list = {im_list.name}';
out_dir = 'Z:\HEproject\data\normalization_512_jpg_gray\';

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
parfor i = 1:length(im_list)
    im_name = im_list{i};
    RGB = imread(fullfile(im_dir,im_name));
    I = rgb2gray(RGB);
    imwrite(I,fullfile(out_dir,im_name));
end
