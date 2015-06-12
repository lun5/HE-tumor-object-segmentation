imdir = 'Z:\Lung\tiles';
imnames = {'290677_im29_stx21583_sty2002.jpg','13fsu-8_im30_stx10361_sty4003.jpg',...
    '13fsu-8_im56_stx14505_sty8005.jpg','13fsu-8_im79_stx12433_sty12007.jpg',...
    '13fsu-8_im101_stx8289_sty16009.jpg','13fsu-2_im37_stx4061_sty4003.jpg',...
    '13fsu-2_im55_stx6091_sty6004.jpg','13fsu-18_im71_stx20351_sty8005.jpg',...
    '13fsu-18_im84_stx16281_sty10006.jpg','13fsu-18_im155_stx8141_sty20011.jpg',...
    '13fsu-18_im131_stx20351_sty16009.jpg'};

numImages = length(imnames);

for i = 1:numImages
    imname = fullfile(imdir,imnames{i});
    raw_image = imread(imname);
    ndown = 4;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
    %figure; imshow(raw_image);rect = getrect;dz_im = imcrop(raw_image,rect);
    I = double(dz_im);figure; imshow(I/255);size(I)
    
    opts_affinity = setEnvironment_affinity;
    [A,im_sizes] = getW(I,opts_affinity);
    opts_clustering = setEnvironment_clustering;
    [segmented_image, E, E_oriented] = graphSegmentation(A,im_sizes,I,opts_clustering);
end

    