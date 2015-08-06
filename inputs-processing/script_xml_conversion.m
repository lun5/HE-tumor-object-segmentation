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

xmlfilename = '2ALe5NgRyfnpo.xml';
folder_name = 'best_images_july30';
new_xml = xmlLabelMe(fullfile(data_dir,'Annotations_ImageScope',xmlfilename),folder_name,'luong');

fileId = fopen(fullfile(HOMEANNOTATIONS,folder_name,xmlfilename),'w');
fprintf(fileId,'%s',new_xml);
fclose(fileId);

database = LMdatabase(fullfile(HOMEANNOTATIONS,folder_name));
[D,j] = LMquery(database, 'folder',folder_name);
[img, seg, names] = LM2segments(D(1), [], HOMEIMAGES, HOMELMSEGMENTS);