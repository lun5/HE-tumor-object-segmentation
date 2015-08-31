% demo converting from annotations to masks
% then masks to boundaries for the evaluation 
clearvars; close all;
addpath(genpath(pwd));
% data_dir = 'Z:\HEproject\data';
% HOMELMSEGMENTS = fullfile(data_dir,'groundTruth_512_512');
% gt_img_dir = 'Z:\TextonInput\all_gt';%fullfile(data_dir,'gt_images_textonBoost');
data_dir = '/Users/lun5/Research/data/';
HOMELMSEGMENTS = fullfile(data_dir,'groundTruth','groundTruth_512_512');%'_reannotated','best_images_july30');
gt_img_dir = fullfile(data_dir,'groundTruth','groundTruth_512_display');
if ~ exist(gt_img_dir,'dir');
    mkdir(gt_img_dir);
end

filenames = dir(fullfile(HOMELMSEGMENTS,'*.mat'));
filenames = {filenames.name}';
numfiles = length(filenames);
components = {'carcinoma','vessel','fat','stroma','lymphocyte','duct',...
    'white space','atypical ductal hyperplasia'};%,'benign terminal lobular unit','unsure'};
num_comps = length(components);
%colors = [1 0 0; 0 1 0; 0 0 1; 1 0.55 0;1 0 1;0 1 1;0 0 0;1 1 0; 0.58 0 0.83;0.72 0.53 0.04];
colors = zeros(num_comps,3);
% mapping colors

for id = 1: length(components)
   colors(id,:) = mapLabel2Color(id-1);
end

%colors = uint8(colors.*255);
parfor i = 1:numfiles
   fname = filenames{i};
   tmp = load(fullfile(HOMELMSEGMENTS,fname));
   gt = tmp.groundTruth{1};
   seg = gt.Segmentation;
   num_objects = length(gt.names);
   for j = 1:num_objects
      obj_id = ismember(components, gt.names{j});
      if sum(obj_id) == 0
          seg(seg == j) = 0;
      else
        seg(seg == j) = find(obj_id) + 20;
      end
   end
   seg = seg - 20;
   seg(seg < 0) = 0;
   seg(seg == 0) = 4;
   gt_img = label2rgb(seg,colors./255);
   imwrite(gt_img,fullfile(gt_img_dir,[fname(1:end-4) '.bmp']));
   fprintf('Done with file %s\n',fname(1:end-4));
end

% filenames = dir(fullfile(HOMELMSEGMENTS,'*.mat'));
% filenames = {filenames.name}';
% numfiles = length(filenames);
% tiles_dir = 'Z:\TilesForLabeling_tiff_renamed';
% new_tiles_dir = 'Z:\TextonInput\all_images';
% parfor i = 1:numfiles
%     fname = filenames{i};
%     im = imread(fullfile(tiles_dir,[fname(1:end-4) '.tif']));
%     im = im(1:4:end,1:4:end,:);
%     %imwrite(im,fullfile(new_tiles_dir,[fname(1:end-4) '.bmp']));
%     imwrite(im,fullfile(new_tiles_dir,[fname(1:end-4) '.tif']),'Compression','none');
%     fprintf('Done with file %s\n',fname(1:end-4));
% end