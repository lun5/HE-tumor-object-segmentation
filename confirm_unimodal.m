addpath(genpath(pwd));
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
% tiles_dir = fullfile(sourcedir,'TilesForLabeling');
% 
% rotation_matrix = load(fullfile(pwd,'rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
% rotation_matrix = rotation_matrix.rotation_matrix;
% % I = imread(fullfile(tiles_dir, 'mws09-778a_12288_12288_2048_2048.tif'));
% I = imread(fullfile(tiles_dir, 'fFwTGXYlhYNa.tif'));
% figure; imshow(I);
% I1 = imcrop; imshow(I1);
% raw_image = double(I1);
% r = raw_image(:,:,1)./255; g = raw_image(:,:,2)./255; b = raw_image(:,:,3)./255;
% rotated_coordinates = rotation_matrix*double([r(:)'; g(:)'; b(:)']);
% theta = atan2(rotated_coordinates(3,:),rotated_coordinates(2,:));
% im_theta = reshape(theta,size(r));

Nsamples = 10000;
opts.sig = 5;
F = sampleF(im_theta,Nsamples,opts);

% figure; imagesc(im_theta); 
% colormap(hsv); colorbar('southoutside'); title('Hue');
% axis equal; axis off; axis tight;
% 
% figure; scatter(F(:,1),F(:,2));
% axis square; axis([-pi pi -pi pi]);
% set(gcf,'color','white') % White background for the figure.
% 
% %figure;[h,xg,yg]=smoothhist2D(F,5,{[-pi:0.1:pi],[-pi:0.1:pi]},.05); 
% figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',2,'points');
% axis equal; axis tight;
% figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1);
% axis equal; axis tight;
figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');

[ params,posterior_probs, prior_probs] = moVM2D(F, 3);
[p,phi,psi] = circ_bvmpdf([],[], 0,0,1,1,1);
i = 1;
[x,y] = meshgrid(phi,psi);

kappa3_hat(:) = 0.5;% = kappa1_hat.*kappa2_hat./(kappa1_hat+kappa2_hat)-2;
figure;
for i = 1:k
p =  circ_bvmpdf(x,y, mu_hat(i), nu_hat(i), kappa1_hat(i),  kappa2_hat(i),kappa3_hat(i));
p = reshape(p,size(x));
contour(x,y,p,'ShowText','on'); axis square; hold on; axis tight;
%figure; surf(x,y,p);axis square; axis tight;
end
set(gcf,'color','white') 
