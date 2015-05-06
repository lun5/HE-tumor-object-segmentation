% addpath(genpath(pwd));
close all; 
clearvars;

% mu = [2.5 2.5 2.5 0 0 -1.5];
% nu = [2.5 0 -1.5 0 -1.5 -1.5];
% kappa1 = [5 5 5 10 10 15];
% kappa2 = [5 10 15 10 15 15];
% kappa3 = [0 0 0 0 0 0];
% n1 = 500; n2 = 100; n3 = 10; n4 = 800; n5 = 150; n6 = 300;
% % 
% alpha1 = circ_bvmrnd([mu(1) nu(1) kappa1(1) kappa2(1) kappa3(1)], n1);
% alpha2 = circ_bvmrnd([mu(2) nu(2) kappa1(2) kappa2(2) kappa3(2)], n2);
% alpha3 = circ_bvmrnd([mu(3) nu(3) kappa1(3) kappa2(3) kappa3(3)], n3);
% alpha4 = circ_bvmrnd([mu(4) nu(4) kappa1(4) kappa2(4) kappa3(4)], n4);
% alpha5 = circ_bvmrnd([mu(5) nu(5) kappa1(5) kappa2(5) kappa3(5)], n5);
% alpha6 = circ_bvmrnd([mu(6) nu(6) kappa1(6) kappa2(6) kappa3(6)], n6);
% 
% alpha = [alpha1; alpha2; alpha3; alpha4; alpha5; alpha6];
% figure; ndhist(alpha(:,1),alpha(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
% axis square;xlabel('phi');ylabel('psi');set(gcf,'color','white');
% % % 
% numData = n1+n2+n3+n4+n5+n6;
% mixtureModel = @(x,y) n1/numData*circ_bvmpdf(x,y,mu(1),nu(1),kappa1(1),kappa2(1),kappa3(1)) + ...
%     n2/numData*circ_bvmpdf(x,y,mu(2),nu(2),kappa1(2),kappa2(2),kappa3(2)) + ...
%     n3/numData*circ_bvmpdf(x,y,mu(3),nu(3),kappa1(3),kappa2(3),kappa3(3)) + ...
%     n4/numData*circ_bvmpdf(x,y,mu(4),nu(4),kappa1(4),kappa2(4),kappa3(4)) + ...
%     n5/numData*circ_bvmpdf(x,y,mu(5),nu(5),kappa1(5),kappa2(5),kappa3(5)) + ...
%     n6/numData*circ_bvmpdf(x,y,mu(6),nu(6),kappa1(6),kappa2(6),kappa3(6));
% [xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
% ppp = mixtureModel(xx,yy);
% ppp = reshape(ppp,size(xx));
% figure;contour(xx,yy,ppp,'ShowText','on');axis square;axis tight;
% set(gcf,'color','white') 
% figure; mesh(xx,yy,ppp);axis square;axis tight;set(gcf,'color','white') 
% [ params,posterior_probs, prior_probs] = mixture_of_bivariate_VM(alpha, 6);
% est_mixtureModel = @(x,y) prior_probs(1)*circ_bvmpdf(x,y,params.mu(1),params.nu(1),params.kappa1(1),params.kappa2(1),params.kappa3(1)) + ...
%     prior_probs(2)*circ_bvmpdf(x,y,params.mu(2),params.nu(2),params.kappa1(2),params.kappa2(2),params.kappa3(2)) + ...
%     prior_probs(3)*circ_bvmpdf(x,y,params.mu(3),params.nu(3),params.kappa1(3),params.kappa2(3),params.kappa3(3)) + ...
%     prior_probs(4)*circ_bvmpdf(x,y,params.mu(4),params.nu(4),params.kappa1(4),params.kappa2(4),params.kappa3(4)) + ...
%     prior_probs(5)*circ_bvmpdf(x,y,params.mu(5),params.nu(5),params.kappa1(5),params.kappa2(5),params.kappa3(5)) + ...
%     prior_probs(6)*circ_bvmpdf(x,y,params.mu(6),params.nu(6),params.kappa1(6),params.kappa2(6),params.kappa3(6));
% [xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
% ppp = est_mixtureModel(xx,yy);
% ppp = reshape(ppp,size(xx));
% figure;contour(xx,yy,ppp,'ShowText','on');axis square;axis tight;set(gcf,'color','white') 
% 
% figure;mesh(xx,yy,ppp);axis square;axis tight;set(gcf,'color','white') 

%%
% Read the image
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
imname = '9uixINHtjjiS.tif';
splitStr = regexp(imname,'\.','split');
raw_image = imread(fullfile(tiles_dir, imname));
I = double(raw_image); 
figure; imshow(raw_image);axis tight;
% calculate the feature
opts_affinity = setEnvironment_affinity;
which_features = opts_affinity.features.which_features;
f_maps = getFeatures(double(I),1,which_features,opts_affinity);
f_map_opphue = f_maps{1};
% sample from the image
Nsamples = 10000;
F = sampleF(f_map_opphue,Nsamples,opts_affinity);  
figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
zlim([0 0.05])
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 
% Estimate the mixture model
% precalculate the intital parameters
alldata = f_map_opphue(:);
[ mu_hat_polar,~, kappa_hat,~, prior_probs] = moVM([cos(alldata(:,1)) sin(alldata(:,1))],3);
init_params.theta_hat = mu_hat_polar;
init_params.kappa_hat = kappa_hat;

[ params,posterior_probs, prior_probs] = mixture_of_bivariate_VM(F, 9, init_params);
est_mixtureModel = @(x,y) prior_probs(1)*circ_bvmpdf(x,y,params.mu(1),params.nu(1),params.kappa1(1),params.kappa2(1),params.kappa3(1)) + ...
    prior_probs(2)*circ_bvmpdf(x,y,params.mu(2),params.nu(2),params.kappa1(2),params.kappa2(2),params.kappa3(2)) + ...
    prior_probs(3)*circ_bvmpdf(x,y,params.mu(3),params.nu(3),params.kappa1(3),params.kappa2(3),params.kappa3(3)) + ...
    prior_probs(4)*circ_bvmpdf(x,y,params.mu(4),params.nu(4),params.kappa1(4),params.kappa2(4),params.kappa3(4)) + ...
    prior_probs(5)*circ_bvmpdf(x,y,params.mu(5),params.nu(5),params.kappa1(5),params.kappa2(5),params.kappa3(5)) + ...
    prior_probs(6)*circ_bvmpdf(x,y,params.mu(6),params.nu(6),params.kappa1(6),params.kappa2(6),params.kappa3(6)) + ...
    prior_probs(7)*circ_bvmpdf(x,y,params.mu(7),params.nu(7),params.kappa1(7),params.kappa2(7),params.kappa3(7)) + ...
    prior_probs(8)*circ_bvmpdf(x,y,params.mu(8),params.nu(8),params.kappa1(8),params.kappa2(8),params.kappa3(8)) + ...
    prior_probs(9)*circ_bvmpdf(x,y,params.mu(9),params.nu(9),params.kappa1(9),params.kappa2(9),params.kappa3(9));
[xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
ppp = est_mixtureModel(xx,yy);
ppp = reshape(ppp,size(xx));
figure;contourf(xx,yy,ppp); colorbar; axis square;axis tight;set(gcf,'color','white') 
xlabel('\phi_A'); ylabel('\phi_B'); set(gca,'FontSize',16);

figure;mesh(xx,yy,ppp);axis square;axis tight;set(gcf,'color','white') 
xlabel('\phi_A'); ylabel('\phi_B'); set(gca,'FontSize',16);
% Q: how do you go from pdf to probability?

mixture_params.params = params;
mixture_params.prior_probs = prior_probs;
[pmi,pJoint,pProd] = evalPMI_theta(F,mixture_params,opts_affinity);





