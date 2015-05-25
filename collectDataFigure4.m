%% Sample the pixel pairs from all images 
tiles_dir = fullfile('Z:\','TilesForLabeling_tiff_renamed'); %window
gt_dir = 'Z:\HEproject\data\groundTruth_2048_2048'; %window

fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths);% 232
opts_affinity = setEnvironment_affinity;
Nsamples = 1000;
parfor j = 1:numImages
    imname = imagepaths{j}; %imname = 'h1402uhfkz.tif';
    im_splitStr = regexp(imname,'\.','split');
    raw_image = imread(fullfile(tiles_dir,imname));
    I = double(raw_image);    
    % sampleF of the RGB image
    [F,p1,p2] = sampleF(f_maps,Nsamples,opts,mask)
    % read the ground truth
    data = load(fullfile(gt_dir,[im_splitStr{1} '.mat']));
    segment = data{1.1}.Segmentation;
    
    % check if the pixels are in the same segments
    indA = sub2ind(size(raw_image),rowSubA,colSubA);
    labelsA = segment(indA);
    indB = sub2ind(size(raw_image),rowSubB,colSubB);
    labelsB = segment(indB);
    labels = (labelsA == labelsB);
    savedata = [indA indB featureA featureB labels];
    fname = fullfile('Z:\HEproject\data\pixelsPair\',[im_splitStr{1} '_pairs.mat']);
    parsave(fname,savedata);
end

