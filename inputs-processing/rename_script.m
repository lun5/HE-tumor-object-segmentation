% %% This script generate file names for the jpeg files
% Luong Nguyen
% 3/17/2015
addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_jpg');
tiles_renamed_dir = fullfile(sourcedir,'TilesForLabeling_jpg_renamed');

if ~exist(tiles_renamed_dir,'dir')
    mkdir(tiles_renamed_dir);
    fileattrib(tiles_renamed_dir,'+w');
end

symbols = ['a':'z' 'A':'Z' '0':'9'];
MAX_ST_LENGTH = 15;MIN_ST_LENGTH = 10;
fileNames = dir(fullfile(tiles_dir,'*.jpg'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths); 

renamed_images = imagepaths;
st = '';
for i = 1:numImages
    stLength = randi([MIN_ST_LENGTH,MAX_ST_LENGTH]);
    nums = randi(numel(symbols),[1 stLength]);
    st = symbols(nums);
    while ismember(st,renamed_images)
        nums = randi(numel(symbols),[1 stLength]);
        st = symbols(nums);
    end    
    renamed_images{i} = st;
    copyfile(fullfile(tiles_dir,imagepaths{i}),fullfile(tiles_renamed_dir,[st '.jpg']));
end
mapping_names = [imagepaths,renamed_images];
T = cell2table(mapping_names,'VariableNames',{'Old_Names','New_Names'});
writetable(T,fullfile(sourcedir,'mappingNames.txt'));

%T = readtable(fullfile(sourcedir,'mappingNames.txt'));
C = table2cell(T);
renamed_images = C(:,2);

fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths); 

for i = 1:numImages
    copyfile(fullfile(tiles_dir,imagepaths{i}),...
        fullfile(tiles_renamed_dir,[renamed_images{i} '.tif']));
end

