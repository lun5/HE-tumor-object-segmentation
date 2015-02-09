addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling');
rotation_matrix = load(fullfile(pwd,'rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
rotation_matrix = rotation_matrix.rotation_matrix;
I = imread(fullfile(tiles_dir, 'mws09-778a_12288_12288_2048_2048.tif'));
figure; imshow(I);
I1 = imcrop; imshow(I1);
raw_image = double(I1);
r = raw_image(:,:,1)./255; g = raw_image(:,:,2)./255; b = raw_image(:,:,3)./255;
rotated_coordinates = rotation_matrix*double([r(:)'; g(:)'; b(:)']);
theta = atan2(rotated_coordinates(3,:),rotated_coordinates(2,:));
im_theta = reshape(theta,size(r));

Nsamples = 10000;
opts.sig = 3;
F = sampleF(im_theta,Nsamples,opts);

figure; imagesc(im_theta); 
colormap(hsv); colorbar('southoutside'); title('Hue');
axis equal; axis off; axis tight;

figure; scatter(F(:,1),F(:,2));
axis square; axis([-pi pi -pi pi]);
set(gcf,'color','white') % White background for the figure.

% X = [mvnrnd([0 5], [3 0; 0 3], 2000);
%            mvnrnd([0 8], [1 0; 0 5], 2000);
%            mvnrnd([3 5], [5 0; 0 1], 2000)];
% smoothhist2D(X,5,[100, 100],.05);
% smoothhist2D(X,5,[100, 100],[],'surf');
figure;[h,xg,yg]=smoothhist2D(F,5,{[-pi:0.1:pi],[-pi:0.1:pi]},.05); 