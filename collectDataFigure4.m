%% Sample the pixel pairs from all images 
tiles_dir = fullfile('Z:\','TilesForLabeling_tiff_renamed'); %window
gt_dir = 'Z:\HEproject\data\groundTruth_2048_2048\users\lun5\renamed_images'; %window

fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths);% 232
opts_affinity = setEnvironment_affinity;
Nsamples = 1000;
%sample = zeros(Nsamples*numImages,9);
sample = cell(1,numImages);
parfor j = 1:numImages
    imname = imagepaths{j}; %imname = 'h1402uhfkz.tif';
    im_splitStr = regexp(imname,'\.','split');
    raw_image = imread(fullfile(tiles_dir,imname));
    I = double(raw_image);    
    % sampleF of the RGB image
    [F,p1,p2] = sampleF(I,Nsamples,opts_affinity);
    % read the ground truth
    data = load(fullfile(gt_dir,[im_splitStr{1} '.mat']));
    segment = data.groundTruth{1,1}.Segmentation;
    im_size = size(raw_image(:,:,1));
    % check if the pixels are in the same segments
    ii = sub2ind(im_size,p1(:,1),p1(:,2));
    jj = sub2ind(im_size,p2(:,1),p2(:,2));
    labelsA = segment(ii);
    labelsB = segment(jj);
    labels = (labelsA == labelsB);
    sample{j} = [ii jj F labels];
    %sample(((j-1)*Nsamples+1):((j-1)*Nsamples + nAcceptedData),:) = [ii jj F labels];
    display(['finish with image ', imname]);
end

fname = fullfile('Z:\HEproject\data\pixelPairs\','sample_pairs_labels.mat');
save(fname,'sample');
numPos = 0; numNeg = 0; total = 0;
for j = 1:numImages
    data = sample{j};
    total = total + size(data,1);
    numPos = numPos + sum(data(:,end));
    numNeg = numNeg + sum(data(:,end) == 0);
end

%pos: 224368; neg: 6135; tot: 230530 --> so much more positive than
%negative



