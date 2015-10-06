%% convert TIFF to JPG lossless matlab
% Luong Nguyen
% 10/5/15

tiff_dir = 'Z:\Tiles_512';
jpg_dir = 'Z:\Tiles_512_jpg';

if ~exist(jpg_dir,'dir')
    mkdir(jpg_dir)
end

im_list = dir(fullfile(tiff_dir,'*.tif'));
im_list = {im_list.name}';
num_images = length(im_list);

for i =1:num_images
   im_name = im_list{i}(1:end-4);
   I = imread(fullfile(tiff_dir,im_list{i}));
   imwrite(I,fullfile(jpg_dir,[im_name '.jpg']));
   fprintf('Done with image %s\n',im_name);
end