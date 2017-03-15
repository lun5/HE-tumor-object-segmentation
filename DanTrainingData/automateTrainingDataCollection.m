% script to extract training data 
% Luong Nguyen 9/14/15

%svs_name = 'tp10-867-1';
svs_name = 'TCGA-05-4244-01A-01-BS1';
%im_dir = 'Z:\TilesForLabeling';
im_dir ='D:\Documents\Tiles_Norm\lung_data\sample_WSI';
%fileList = dir(fullfile(im_dir,[svs_name '*.tif']));
fileList = dir(fullfile(im_dir,[svs_name '_files'], '0','good_tiles','*.jpeg'));
im_names = {fileList.name}';

Hpixels = [];
Epixels = []; 

for i = 1:length(im_names)
    imname = im_names{i};
    I = imread(fullfile(im_dir,[svs_name '_files'], '0','good_tiles',imname));
    [~, H, E] =  extractStainPixels(I,[],[],[1 1.5]);
    Hpixels = cat(1,Hpixels, H);
    Epixels = cat(1,Epixels, E);
end

display('Done');

indx = 1:100:size(Hpixels,1);
figure; scatter3(Hpixels(indx,1),Hpixels(indx,2), Hpixels(indx,3),30,Hpixels(indx,:)./255,'filled');
hold on;scatter3(Epixels(indx,1),Epixels(indx,2), Epixels(indx,3),30,Epixels(indx,:)./255,'filled');