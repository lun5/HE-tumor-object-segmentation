% script to convert xml file from ImageScope to xml for LabelMe
% Luong Nguyen
% 7.30.2015

% List of files
%
LABELME = 'C:\Users\luong_nguyen\Documents\GitHub\LabelMeToolbox';
addpath(genpath(LABELME));
data_dir = 'Z:\TilesForLabeling_bestImages';
HOMEIMAGES = fullfile(data_dir,'Images'); % you can set here your default folder
HOMEANNOTATIONS = fullfile(data_dir,'Annotations'); % you can set here your default folder
HOMELMSEGMENTS = fullfile(data_dir,'groundTruth');

xmlfilename = 'jbaKL4TsEqT.xml';

filenames = dir(fullfile(data_dir,'Annotations_ImageScope','*.xml'));
filenames = {filenames.name}';
folder_name = 'best_images_july30';
parfor i = 1:length(filenames)
    xmlfilename = filenames{i};
    new_xml = xmlLabelMe(fullfile(data_dir,'Annotations_ImageScope',xmlfilename),folder_name,'luong');
    fileId = fopen(fullfile(HOMEANNOTATIONS,folder_name,xmlfilename),'w');
    fprintf(fileId,'%s',new_xml);
    fclose(fileId);
end

database = LMdatabase(fullfile(HOMEANNOTATIONS,folder_name));
[D,j] = LMquery(database, 'folder',folder_name);
[img, seg, names] = LM2segments(D(11), [], HOMEIMAGES, HOMELMSEGMENTS);