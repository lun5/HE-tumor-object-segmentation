addpath(genpath(pwd));
close all; clearvars;
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%tiles_dir = fullfile(sourcedir,'TilesForLabeling');

% tiles_dir = '/Users/lun5/Box Sync/TilesForLabeling';
rotation_matrix = load(fullfile(pwd,'rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
rotation_matrix = rotation_matrix.rotation_matrix;
% I = imread(fullfile(tiles_dir, 'mws09-778a_12288_12288_2048_2048.tif'));
% I = imread(fullfile(tiles_dir, 'fFwTGXYlhYNa.tif'));
I = imread(fullfile(tiles_dir, '2ALe5NgRyfnpo.tif'));
figure; imshow(I);
%I1 = imcrop; imshow(I1);
%raw_image = double(I1);
raw_image = double(I);
r = raw_image(:,:,1)./255; g = raw_image(:,:,2)./255; b = raw_image(:,:,3)./255;
rotated_coordinates = rotation_matrix*double([r(:)'; g(:)'; b(:)']);
theta = atan2(rotated_coordinates(3,:),rotated_coordinates(2,:));
im_theta = reshape(theta,size(r));

Nsamples = 10000;
opts.sig = 7;
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
xlabel('phi'); ylabel('psi'); axis square;set(gca,'FontSize',16);set(gcf,'color','white') 
[ params,posterior_probs, prior_probs] = mixture_of_bivariate_VM(F, 6);
est_mixtureModel = @(x,y) prior_probs(1)*circ_bvmpdf(x,y,params.mu(1),params.nu(1),params.kappa1(1),params.kappa2(1),params.kappa3(1)) + ...
    prior_probs(2)*circ_bvmpdf(x,y,params.mu(2),params.nu(2),params.kappa1(2),params.kappa2(2),params.kappa3(2)) + ...
    prior_probs(3)*circ_bvmpdf(x,y,params.mu(3),params.nu(3),params.kappa1(3),params.kappa2(3),params.kappa3(3)) + ...
    prior_probs(4)*circ_bvmpdf(x,y,params.mu(4),params.nu(4),params.kappa1(4),params.kappa2(4),params.kappa3(4)) + ...
    prior_probs(5)*circ_bvmpdf(x,y,params.mu(5),params.nu(5),params.kappa1(5),params.kappa2(5),params.kappa3(5)) + ...
    prior_probs(6)*circ_bvmpdf(x,y,params.mu(6),params.nu(6),params.kappa1(6),params.kappa2(6),params.kappa3(6));
[xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
ppp = est_mixtureModel(xx,yy);
ppp = reshape(ppp,size(xx));
for numContours = 10:10:50
figure;contour3(xx,yy,ppp,numContours,'ShowText','off');axis square;axis tight;
set(gcf,'color','white');
xlabel('\phi'); ylabel('\psi');set(gca,'FontSize',16);
end

for cl = 1:6
    myfun = @(x,y) prior_probs(cl)*circ_bvmpdf(x,y,params.mu(cl),params.nu(cl),params.kappa1(cl),params.kappa2(cl),params.kappa3(cl)); 
    ppp = myfun(xx,yy);
    ppp = reshape(ppp,size(xx));
    figure;contour(xx,yy,ppp,'ShowText','on');axis square;axis tight;
    xlabel('\phi'); ylabel('\psi'); set(gca,'FontSize',16);
    set(gcf,'color','white') 
end
