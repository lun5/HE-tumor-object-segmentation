% list files with all the k_3
close all; clearvars;
addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
superpixel_input_dir = fullfile(sourcedir,'super_pixels_1');

fileNames = dir(fullfile(superpixel_input_dir,'*_k3'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths);% 420
cd(superpixel_input_dir)
%image_indx = 2:20:numImages;
tic;
parfor j = 1:5%numImages
    imname = imagepaths{j}; 
    im_splitStr = regexp(imname,'\_','split');
    javaCommand = ['java -Xmx2000m -jar ExtractObjectInformation.jar ' ...
        im_splitStr{1} ' 3 5 10 1 128 0'];%        im_splitStr{1} ' 3 3 3 1 128'];
    %tic;
    [status,cmdout] = dos(javaCommand);
    %time_lapsed = toc;
    %sprintf('Finish with file %s in %d seconds\n',im_splitStr{1}, time_lapsed)
    sprintf('Finish with file %s\n',im_splitStr{1})
end
toc