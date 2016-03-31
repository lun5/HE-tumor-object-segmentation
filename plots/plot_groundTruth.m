% %% generate Figure 1
% DATA_DIR = '/Users/lun5/Research/data';
% IMG_DIR = fullfile(DATA_DIR,'Tiles_512');
% GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','well_defined');
% %imname = 'lszomrlgsc5na4q';
% imname = '95f7k8loesyevi';
% img_rgb = imread(fullfile(IMG_DIR,[imname '.tif']));
% img_rgb = double(img_rgb)./255;
% r = img_rgb(:,:,1); g = img_rgb(:,:,2); b = img_rgb(:,:,3);
% tmp = load(fullfile(GT_DIR,[imname '.mat']));
% gt = tmp.groundTruth{1};
% seg = gt.Segmentation;
% names = gt.names;
% nrows = size(seg,1); ncols = size(seg,2);
% bdry = zeros(nrows, ncols);
% num_comps = 0; se = strel('disk',2,4);
% num_classes = length(names);
% for cl = 1:num_classes
%     obj_name = names{cl};
%     if strcmp(obj_name, 'stroma') || strcmp(obj_name, 'white')
%         continue;
%     end
%     CC = bwconncomp(seg==cl);
%     for cp = 1:length(CC.PixelIdxList)
%         % recreate image crop from the mask
%         BW = zeros(CC.ImageSize);
%         BW(CC.PixelIdxList{cp}) = 1;
%         area_mask = sum(BW(:));
%         if area_mask < 1000
%             continue;
%         end
%         num_comps = num_comps + 1;
%         bdry_comp = seg2bdry(BW,'imageSize');
%         bdry_comp = logical(imdilate(bdry_comp,se));
%         bdry(bdry_comp) = num_comps;
%     end
% end
% 
% colors = distinguishable_colors(num_comps);
% for i = 1:num_comps
%     r(bdry == i) = colors(i,1);
%     g(bdry == i) = colors(i,2);
%     b(bdry == i) = colors(i,3);
% end
% 
% output = cat(3,r,g,b);
% imwrite(output,[imname '.tif'],'Resolution',300);

%% generate figure 2
clearvars; close all;
tiles_dir = '/Users/lun5/Research/data/TilesForLabeling_tiff_renamed';
imname = 'lszomrlgsc5na4q';
imname = lower(imname);
raw_image = imread(fullfile(tiles_dir, [imname '.tif']));%figure;imshow(raw_image);
ndown = 1;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
im_rgb = double(dz_im);figure; imshow(im_rgb./255);
nrows = size(im_rgb,1); ncols = size(im_rgb,2);
% angle space
rgb_coords = reshape(im_rgb,[nrows*ncols,3])';
rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
rotation_matrix = rotation_matrix.rotation_matrix;
mask_white = isolateWhite(im_rgb);
indx_white = mask_white(:);
mask_red = isolateRed(im_rgb);
indx_red = mask_red(:);
mu_white = 2.24; kappa_white = 30; % vM mean and concentration of white
rotated_coordinates = rotation_matrix*rgb_coords;
theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
%theta(indx_white) = circ_vmrnd(mu_white, kappa_white, sum(indx_white));
im_theta = reshape(theta,[nrows ncols]);
%figure; imagesc(im_theta); axis off; axis square; colorbar; 
h1 = figure; imagesc(im_theta);
cmap = colormap(hsv); colorbar('southoutside');
axis equal; axis off; axis tight;set(gca,'FontSize',20);
c1 = cmap(1:11,:); c2 = cmap(12:22,:); c3 = cmap(23:32,:);
c4 = cmap(33:43,:); c5 = cmap(44:54,:); c6 = cmap(55:end,:);
cmap_new = [c1;c2;c4;c5;c3;c6];
cut_indx = 11;
cmap_new = [cmap((cut_indx+1):end,:);cmap(1:cut_indx,:)];
h2 = figure; imagesc(im_theta);axis equal; axis off; axis tight;%title('Hue');
colormap(h2,cmap_new);colorbar('southoutside');set(gcf,'color','white');
set(gca,'FontSize',30);
%close(h1);
caxis([-pi pi]); 
print([imname '_theta'],'-dtiff','-r300');

% SIC space
mu_s = 2; sigma_s = 3;
[ sic_coords ] = rgb2sic( rgb_coords, mu_s, sigma_s,[]);% rotation_matrix);
sic1_im = reshape(sic_coords(1,:),[size(im_rgb,1), size(im_rgb,2)]);
sic2_im = reshape(sic_coords(2,:),[size(im_rgb,1), size(im_rgb,2)]);
figure; imagesc(sic1_im); axis off; axis square; colorbar('southoutside');
caxis([-1 1]);set(gca,'FontSize',20);
print([imname '_SIC1'],'-dtiff','-r300');
figure; imagesc(sic2_im); axis off; axis square; colorbar('southoutside');
caxis([-1 1]);set(gca,'FontSize',20);
print([imname '_SIC_2'],'-dtiff','-r300');


