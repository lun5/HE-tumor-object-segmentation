%% Sample the pixel pairs from all images 
tiles_dir = fullfile('Z:\','TilesForLabeling_tiff_renamed'); %window
gt_dir = 'Z:\HEproject\data\groundTruth_2048_2048\users\lun5\renamed_images'; %window

fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths);% 232
opts_affinity = setEnvironment_affinity;
Nsamples = 20000;
%sample = zeros(Nsamples*numImages,9);
sample = cell(1,numImages);
parfor j = 1:numImages
    imname = imagepaths{j}; %imname = 'h1402uhfkz.tif';
    im_splitStr = regexp(imname,'\.','split');
    raw_image = imread(fullfile(tiles_dir,imname));
    I = double(raw_image);    
    % sampleF of the RGB image
    [F,ii,jj] = sampleF(I,Nsamples,opts_affinity);
    % read the ground truth
    data = load(fullfile(gt_dir,[im_splitStr{1} '.mat']));
    segment = data.groundTruth{1,1}.Segmentation;
    im_size = size(raw_image(:,:,1));
    % check if the pixels are in the same segments
    labelsA = segment(ii);
    labelsB = segment(jj);
    labels = (labelsA == labelsB);
    sample{j} = [ii jj F labels];   
    display(['finish with image ', imname]);
end

fname = fullfile('Z:\HEproject\data\pixelPairs\','RGBsample_pairs_labels.mat');
save(fname,'sample');
numPos = 0; numNeg = 0; 
for j = 1:numImages
    data = sample{j};
    numPos = numPos + sum(data(:,end));
    numNeg = numNeg + sum(data(:,end) == 0);
end

%numPos  4476509 ~ 97.3%
%numNeg  123930 ~ 2.7%
%total  4600439
%pos: 224368; neg: 6135; tot: 230530 --> so much more positive than
%negative

outputDir = 'Z:\HEproject\data\pixelPairs\';
%which_features = 'luminance';
%which_features = 'hue opp';
which_features = 'color';
inputFile = fullfile(outputDir,'RGBsample_pairs_labels.mat');
processInputClassifier(inputFile, outputDir, which_features);



