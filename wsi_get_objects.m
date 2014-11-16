% gather training data for segmentation
function wsi_get_objects(svs_fname, datadir, resultdir) 

fileNames = dir(fullfile(datadir,[svs_fname '*.' 'tif']));
imagepaths = {fileNames.name}';

numImages = length(imagepaths);% 420
numImagesQualified = 0;

for i = 1:numImages
   
    imname = imagepaths{i}; 
    raw_image = imread(fullfile(datadir,imname));
    imshow(raw_image);
    choice = menu('Image qualified for traning?','white','qualified');
    while choice == 2 
        rect = getrect;
        im = imcrop(raw_image, rect);
        imwrite(im,fullfile(resultdir,[svs_fname, 'gland', num2str(numImagesQualified),'.tif']));
        numImagesQualified = numImagesQualified + 1;                                
        choice = menu('Image qualified for traning?','white','qualified');
    end

    % only consider up to 20 qualified images
    if numImagesQualified > 40
        break;
    end        
end

end