% script to collect data tiles
sourcedir = 'D:\Documents\Tiles_Norm\lung_data\sample_WSI';
svs_fnames = dir(fullfile(sourcedir,'*.tif'));
svs_fnames = {svs_fnames.name}';
num_svs = length(svs_fnames);
% folder storing the tiles of interest
training_dir = fullfile(sourcedir,'ColorsTrainingData');

if ~exist(training_dir,'dir');mkdir(training_dir);end
%numtiles = 10;

for i = 1:num_svs
    svs_fname = svs_fnames{i}(1:end-4);   
    tiles_dir = fullfile(sourcedir,[svs_fname '_files'],'0','good_tiles');
    fileNames = dir(fullfile(tiles_dir,'*.jpeg'));
    imagepaths = {fileNames.name}';

    numImages = length(imagepaths);% 420
  
    %training_data = zeros(3,numImages*nrepeats*2);
    training_data_purple = zeros(3, 1e6);
    training_data_pink = zeros(3, 1e6);

    count_purple = 1; % number of training examples
    count_pink = 1;
    for j = 1:numImages
        imname = imagepaths{j}; 
        raw_image = imread(fullfile(tiles_dir,imname));

        % double the size of training_data_purple if needed
        if count_purple + 1e3 > size(training_data_purple,2)
            temp = training_data_purple;
            training_data_purple = zeros(3,size(training_data_purple,2)*2);
            training_data_purple(:, size(temp,2)) = temp;
        end
    
        % double the size of training_data_pink if needed
        if count_pink + 1e3 > size(training_data_pink,2)
            temp = training_data_pink;
            training_data_pink = zeros(3,size(training_data_pink,2)*2);
            training_data_pink(:, size(temp,2)) = temp;
        end
         
        % collect purple data
        figure; imshow(raw_image); 
        zoom on;pause(); zoom off; %press enter to get out of the zoom mode
        rect = getrect;zoom out;
        I_purple = imcrop(raw_image,rect);        
        %I_purple = imcrop(raw_image);
        npixels = numel(I_purple(:,:,1));
        purple_image_crop_rgb = [reshape(I_purple(:,:,1),1,npixels);...
            reshape(I_purple(:,:,2),1,npixels);reshape(I_purple(:,:,3),1,npixels)];
        training_data_purple(:,count_purple:count_purple+npixels-1) = purple_image_crop_rgb;
        count_purple = count_purple + npixels;
    

        % collect pink data
        zoom on;pause(); zoom off; %press enter to get out of the zoom mode
        rect = getrect;zoom out;
        I_pink = imcrop(raw_image,rect);
        npixels = numel(I_pink(:,:,1));
        pink_image_crop_rgb = [reshape(I_pink(:,:,1),1,npixels);...
            reshape(I_pink(:,:,2),1,npixels);reshape(I_pink(:,:,3),1,npixels)];
        training_data_pink(:,count_pink:count_pink+npixels-1) = pink_image_crop_rgb;
        count_pink = count_pink + npixels;
        % select pink
        close all;
    end
    % remove excessive space in the pink/purple arrays
    training_data_purple(:, count_purple+1:end) = [];
    training_data_pink(:, count_pink+1:end) = [];
    % save the training data for that image
    save([training_dir filesep splitStr{1} 'training_purple.mat'],'training_data_purple');
    save([training_dir filesep splitStr{1} 'training_pink.mat'],'training_data_pink');
end
