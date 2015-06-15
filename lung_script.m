%imdir = 'Z:\Lung\tiles';
imdir = '/Users/lun5/Research/data/Lung_tiles';
imnames = {'290677_im29_stx21583_sty2002.jpg','13fsu-8_im30_stx10361_sty4003.jpg',...
    '13fsu-8_im56_stx14505_sty8005.jpg','13fsu-8_im79_stx12433_sty12007.jpg',...
    '13fsu-8_im101_stx8289_sty16009.jpg','13fsu-2_im37_stx4061_sty4003.jpg',...
    '13fsu-2_im55_stx6091_sty6004.jpg','13fsu-18_im71_stx20351_sty8005.jpg',...
    '13fsu-18_im84_stx16281_sty10006.jpg','13fsu-18_im155_stx8141_sty20011.jpg',...
    '13fsu-18_im131_stx20351_sty16009.jpg'};
imnames = {'13fsu-18_im71_stx20351_sty8005.jpg','13fsu-18_im84_stx16281_sty10006.jpg'};

result_dir = fullfile(imdir,'results');
if ~exist(result_dir,'dir')
    mkdir(result_dir)
end

numImages = length(imnames);
opts_affinity = setEnvironment_affinity;
which_features = opts_affinity.features.which_features{1};
feature_result_dir = fullfile(result_dir,which_features);
if ~exist(feature_result_dir,'dir')
    mkdir(feature_result_dir)
end

opts_clustering = setEnvironment_clustering;
for i = 1:numImages
    %tic;
    imname = fullfile(imdir,imnames{i});
    raw_image = imread(imname);
    ndown = 4;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
    %figure; imshow(raw_image);rect = getrect;dz_im = imcrop(raw_image,rect);
    I = double(dz_im);%  figure; imshow(I/255);size(I)
    [A,im_sizes] = getW(I,opts_affinity);
    [segmented_image, E, E_oriented] = graphSegmentation(A,im_sizes,I,opts_clustering);
    imwrite(mat2gray(1-E),fullfile(feature_result_dir,[imnames{i}(1:end-4) '_edgemap.tiff']),'Resolution',300);
    imwrite(segmented_image./255,fullfile(feature_result_dir,[imnames{i}(1:end-4) '_segmented.tiff']),'Resolution',300); 
    %close all;
    %elapsed_time = toc; 
    fprintf('done with image %s\n',imnames{i}(1:end-4));
end

    