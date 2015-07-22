%% Create figure to illustrate the ROI ranking idea
% Luong Nguyen 07/20/2015

%% Use image tp09-1313-1 and tp09-1313-2
data_dir = 'Z:\HEproject\data';
%wsi_images = {'tp09-1313-1','tp09-1313-2'};
file_mapping = readtable('Z:\mappingNames.txt');
old_names = file_mapping.Old_Names;
new_names = file_mapping.New_Names;
split_old_names = cellfun(@(x) regexp(x,'_','split'), old_names,'UniformOutput',false);
split_old_names = cat(1,split_old_names{:});
wsi_images = unique(split_old_names(:,1));

%% extract the components
%% query the images
%clearvars; close all;
LABELME = 'C:\Users\luong_nguyen\Documents\GitHub\LabelMeToolbox';
addpath(genpath(LABELME));
HOMEIMAGES = fullfile(data_dir,'Images'); % you can set here your default folder
HOMEANNOTATIONS = fullfile(data_dir,'Annotations'); % you can set here your default folder
HOMELMSEGMENTS = fullfile(data_dir,'groundTruth_fake');
mixture1d_dir = 'Z:\mixture_von_mises\same_rot_renamed_images_FreezeAll';
tiles_dir = 'Z:\TilesForLabeling_tiff_renamed';

folderlist = {'renamed_images'};
database = LMdatabase(fullfile(HOMEANNOTATIONS,'users','lun5','renamed_images'));

for i = 1:length(wsi_images)
   T = tic;
   output_dir = fullfile(data_dir,['tissue_components_' wsi_images{i}]);
   if ~exist(output_dir,'dir')
       mkdir(output_dir);
   end
   [ind_wsi,~] = ismember(split_old_names(:,1),wsi_images{i});
   new_names_wsi = new_names(ind_wsi);
   %% extract the component
   num_images = length(new_names_wsi);
   parfor j = 1:num_images
       imname = new_names_wsi{j}; imname = lower(imname);
       [dd,~] = LMquery(database, 'filename',[imname '.jpg']);
       [~, seg, names] = LM2segments(dd, [], HOMEIMAGES, HOMELMSEGMENTS);
       num_classes = length(names);
       img = double(imread(fullfile(tiles_dir,[imname '.tif'])))./255;
       img_brighten = imadjust(img,[0 0 0;1 1 1],[],0.1);
       for cl = 1:num_classes
           obj_name = names{cl};
           CC = bwconncomp(seg==cl);
           for cp = 1:length(CC.PixelIdxList)
               % recreate image crop from the mask
               BW = zeros(CC.ImageSize);
               BW(CC.PixelIdxList{cp}) = 1;
               bmap = logical(seg2bdry(BW,'imageSize'));
               se = strel('ball',3,3); 
               bdry = imdilate(uint8(bmap),se);
               BW = repmat(BW,[1 1 3]);
               img_combined = (img.*(BW) + img_brighten.*(~BW)).*repmat(~bdry,[1 1 3]);
               img_combined = img_combined + double(cat(3,bdry,zeros(size(bdry,1),size(bdry,2),2)));
               imgCombined_name = fullfile(output_dir,[imname '_obj' num2str(cp)...
                   '_' obj_name '.tif']);
               imwrite(img_combined,imgCombined_name,'Resolution',300);          
           end
       end
   end
   t = toc(T);
   fprintf('done with image %s in %d seconds\n',wsi_images{i},t);
end