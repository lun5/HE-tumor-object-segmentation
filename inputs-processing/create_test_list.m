%tiles_dir = 'Z:\HEproject\normalization_512';
%tiles_dir = 'Z:\HEproject\normalization_2048';
tiles_dir = 'Z:\HEproject\data\GroundTruth\groundTruth_512_512';
%tiles_dir_test = 'Z:\HEproject\Test';
tiles_dir_test ='Z:\TilesForLabeling_bestImages\groundTruth_512_512_reannotated_Oct\best_images_july30';

filelist = dir(fullfile(tiles_dir_test,'*.mat'));
filelist = {filelist.name}';

if ~exist(fullfile(tiles_dir,'well_defined'),'dir')
    mkdir(fullfile(tiles_dir,'well_defined'));
end

if ~exist(fullfile(tiles_dir,'invasive'),'dir')
    mkdir(fullfile(tiles_dir,'invasive'));
end

file_ext = '.mat';
for i = 1:length(filelist)
    imname = filelist{i}(1:end-4);
    imname = lower(imname);
    if exist(fullfile(tiles_dir,'invasive',[imname file_ext]),'file')
        movefile(fullfile(tiles_dir,'invasive',[imname file_ext]), fullfile(tiles_dir,'well_defined',[imname file_ext]));
    end
    %movefile(fullfile(tiles_dir,[imname '.mat']), fullfile(tiles_dir,'well_defined',[imname '.mat']));
end