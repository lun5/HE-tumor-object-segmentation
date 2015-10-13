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

%xmlfilename = 'jbaKL4TsEqT.xml';
filenames = dir(fullfile(data_dir,'Annotations_ImageScope','*.xml'));
filenames = {filenames.name}';
folder_name = 'best_images_july30';

parfor i = 1:length(filenames)
    xmlfilename = filenames{i};
    new_xmlfilename = fullfile(HOMEANNOTATIONS,folder_name,xmlfilename);
    if ~exist(new_xmlfilename,'file')
        T = tic;
        new_xml = xmlLabelMe(fullfile(data_dir,'Annotations_ImageScope',xmlfilename),folder_name,'luong');
        fileId = fopen(fullfile(HOMEANNOTATIONS,folder_name,xmlfilename),'w');
        fprintf(fileId,'%s',new_xml);
        fclose(fileId);
        elapsed_time = toc(T); fprintf('\n Done with %s in %.1f seconds\n',xmlfilename, elapsed_time);
    end
end

database = LMdatabase(fullfile(HOMEANNOTATIONS,folder_name));
[D,j] = LMquery(database, 'folder',folder_name);

bdry_im_dir = fullfile(data_dir,'bdry_im');
seg_im_dir = fullfile(data_dir,'segmentation_display');
if ~exist(bdry_im_dir,'dir')
    mkdir(bdry_im_dir);
end

if ~exist(seg_im_dir,'dir')
    mkdir(seg_im_dir);
end
Nimages = length(D);
parfor ndx = 1:Nimages
    if ~isfield(D(ndx).annotation, 'object')
        continue; 
    end
    im_name = D(ndx).annotation.filename(1:end-4);
    fileseg = fullfile(HOMELMSEGMENTS, D(ndx).annotation.folder, [im_name '.mat']);
    if ~exist(fileseg,'file')
        T = tic;
        [img, seg, names] = LM2segments(D(ndx), [], HOMEIMAGES, HOMELMSEGMENTS);
        mult = 4;
        seg = seg(1:mult:end,1:mult:end);img = imresize(img,size(seg));
        %[annotation, img] = LMread(D,ndx,HOMEIMAGES);
        % plot the annotations %figure;LMplot(annotation, img)
        bmap = logical(seg2bdry(seg,'imageSize'));  %figure; imshow(mat2gray(bmap));
        groundTruth = cell(1,1);
        groundTruth{1,1} = struct('Segmentation',seg,'Boundaries',bmap);
        groundTruth{1,1}.names = names;
        parsave_custom(fileseg,groundTruth,'groundTruth');
        % change thickness of edges
        edge_map = imdilate(bmap, strel('disk',1));
        edge_map_im = img.*uint8(repmat(~edge_map,[1 1 3]));
        imwrite(uint8(edge_map_im),fullfile(bdry_im_dir,[im_name, '_bdry.tif']),'Resolution',300);
        imwrite(uint8(label2rgb(seg)),fullfile(seg_im_dir,[im_name, '_seg.tif']),'Resolution',300);
        elapsed_time = toc(T); fprintf('\n Done with %s in %.1f seconds\n',im_name, elapsed_time);
    end
end

disp('Done');