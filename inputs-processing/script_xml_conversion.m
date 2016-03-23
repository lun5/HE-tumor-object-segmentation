% script to convert xml file from ImageScope to xml for LabelMe
% Luong Nguyen
% 7.30.2015

% List of files
%
LABELME = 'C:\Users\luong_nguyen\Documents\GitHub\LabelMeToolbox';
addpath(genpath(LABELME));
data_dir = 'Z:\TilesForLabeling_bestImages';
HOMEIMAGES = fullfile(data_dir,'Images'); % you can set here your default folder
%HOMEIMAGES = 'Z:\TilesForLabeling_jpg_renamed';
HOMEANNOTATIONS = fullfile(data_dir,'Annotations'); % you can set here your default folder
%HOMEANNOTATIONS = 
HOMELMSEGMENTS = fullfile(data_dir,'groundTruth');
if ~exist(fullfile(HOMELMSEGMENTS,'GT_separateComp'),'dir')
    mkdir(fullfile(HOMELMSEGMENTS,'GT_separateComp'))
end

if ~exist(fullfile(HOMELMSEGMENTS,'bdry_im'),'dir')
    mkdir(fullfile(HOMELMSEGMENTS,'bdry_im'))
end
%xmlfilename = 'jbaKL4TsEqT.xml';
filenames = dir(fullfile(data_dir,'Annotations_ImageScope','*.xml'));
filenames = {filenames.name}';
folder_name = 'best_images_july30';
% 
% parfor i = 1:length(filenames)
%     xmlfilename = filenames{i};
%     new_xmlfilename = fullfile(HOMEANNOTATIONS,folder_name,xmlfilename);
%     if ~exist(new_xmlfilename,'file')
%         T = tic;
%         new_xml = xmlLabelMe(fullfile(data_dir,'Annotations_ImageScope',xmlfilename),folder_name,'luong');
%         fileId = fopen(fullfile(HOMEANNOTATIONS,folder_name,xmlfilename),'w');
%         fprintf(fileId,'%s',new_xml);
%         fclose(fileId);
%         elapsed_time = toc(T); fprintf('\n Done with %s in %.1f seconds\n',xmlfilename, elapsed_time);
%     end
% end

%database = LMdatabase(fullfile(HOMEANNOTATIONS,folder_name));
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
imagesize = [2048 2048];
for ndx = 1:Nimages
    im_name = D(ndx).annotation.filename(1:end-4);
    fileseg = fullfile(HOMELMSEGMENTS,'GT_separateComp', [im_name '.mat']);
    imgtmp = zeros(imagesize);
    annotationtmp = LMsortlayers(D(ndx).annotation, imgtmp);
    [S_instances, classes] = LMobjectmask(annotationtmp, size(imgtmp));
    if ~isempty(S_instances) && ~exist(fileseg,'file') 
        % Assign labels taking into account occlusions!
        Mclasses = zeros([imagesize(1) imagesize(2)]);
        for k = size(S_instances,3):-1:1;
            S_instances(:,:,k) = k*S_instances(:,:,k);
            Mclasses = Mclasses+(Mclasses==0).*S_instances(:,:,k);
        end
        seg  = Mclasses.^0.5;
        mult = 4;
        seg = seg(1:mult:end,1:mult:end);img = imresize(img,size(seg));
        %[annotation, img] = LMread(D,ndx,HOMEIMAGES);
        % plot the annotations %figure;LMplot(annotation, img)
        bmap = logical(seg2bdry(seg,'imageSize'));  %figure; imshow(mat2gray(bmap));
        groundTruth = cell(1,1);
        groundTruth{1,1} = struct('Segmentation',seg,'Boundaries',bmap);
        %groundTruth{1,1}.names = names;
        parsave_custom(fileseg,groundTruth,'groundTruth');
        fprintf('Done with image %s\n',im_name);
    end
end
%HOMELMSEGMENTS = fullfile(data_dir,'groundTruth_512_512_reannotated_Oct');
% for ndx = 1:Nimages
% %     if ~isfield(D(ndx).annotation, 'object')
% %         continue; 
% %     end
%     im_name = D(ndx).annotation.filename(1:end-4);
%     fileseg = fullfile(HOMELMSEGMENTS, D(ndx).annotation.folder, [im_name '.mat']);
%     %if ~exist(fileseg,'file')
%         T = tic;
%         [img, seg, names] = LM2segments(D(ndx), [], HOMEIMAGES, HOMELMSEGMENTS);      
%         mult = 4;
%         seg = seg(1:mult:end,1:mult:end);img = imresize(img,size(seg));
%         %[annotation, img] = LMread(D,ndx,HOMEIMAGES);
%         % plot the annotations %figure;LMplot(annotation, img)
%         bmap = logical(seg2bdry(seg,'imageSize'));  %figure; imshow(mat2gray(bmap));
%         groundTruth = cell(1,1);
%         groundTruth{1,1} = struct('Segmentation',seg,'Boundaries',bmap);
%         groundTruth{1,1}.names = names;
%         %parsave_custom(fileseg,groundTruth,'groundTruth');
%         % change thickness of edges
%         edge_map = imdilate(bmap, strel('disk',1));
%         edge_map_im = img.*uint8(repmat(~edge_map,[1 1 3]));
%         imwrite(uint8(edge_map_im),fullfile(bdry_im_dir,[im_name, '_bdry.tif']),'Resolution',300);
%         imwrite(uint8(label2rgb(seg)),fullfile(seg_im_dir,[im_name, '_seg.tif']),'Resolution',300);
%         elapsed_time = toc(T); fprintf('\n Done with %s in %.1f seconds\n',im_name, elapsed_time);
%     %end
% end

%disp('Done');
% gt_dir = 'Z:\HEproject\data\GroundTruth\coarse_fine_GT_512_512\well_defined';
% filenames = dir(fullfile(gt_dir,'*.mat'));
% filenames = {filenames.name}';
% outdir = 'Z:\HEproject\data\GroundTruth\coarse_fine_GT_512_512_new\well_defined';
% 
% if ~exist(outdir,'dir')
%     mkdir(outdir);
% end
% 
% parfor i = 1:length(filenames)
%    imname = filenames{i}(1:end-4);
%    tmp = load(fullfile(gt_dir,[imname '.mat']));
%    groundTruth = tmp.groundTruth;
%    bmap = groundTruth{1,1}.Boundaries;
%    bmap = imdilate(bmap,strel('disk',2));
%    groundTruth{1,1}.Boundaries = bmap;
%    parsave_custom(fullfile(outdir,[imname '.mat']),groundTruth,'groundTruth');
% end
% 
% mkdir(fullfile(outdir),'well_defined')
% mkdir(fullfile(outdir),'invasive');
% 
% gt_dir = 'Z:\HEproject\data\GroundTruth\coarse_fine_GT_512_512\well_defined';
% filenames = dir(fullfile(gt_dir,'*.mat'));
% filenames = {filenames.name}';
% outdir = 'Z:\HEproject\data\GroundTruth\coarse_fine_GT_512_512_new\';
% 
% parfor i = 1:length(filenames)
%    imname = filenames{i}(1:end-4);
%    copyfile(fullfile(outdir,'well_defined',[imname '.mat']), fullfile(outdir,'all_files'))
% end
% 
