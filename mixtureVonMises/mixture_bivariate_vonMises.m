% addpath(genpath(pwd));
close all; 
clearvars;

mu = [0 2];
nu = [0 2];
kappa1 = [5 15];
kappa2 = [15 5];
kappa3 = [0 0];
n1 = 200; n2 = 100; n3 = 300;
% 
alpha1 = circ_bvmrnd([mu(1) nu(1) kappa1(1) kappa2(1) kappa3(1)], n1);
alpha2 = circ_bvmrnd([mu(2) nu(2) kappa1(2) kappa2(2) kappa3(2)], n2);
alpha3 = circ_bvmrnd([mu(2) nu(1) kappa1(2) kappa2(2) kappa3(2)], n3);
alpha = [alpha1; alpha2; alpha3];
figure; ndhist(alpha(:,1),alpha(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
axis square;xlabel('phi');ylabel('psi');set(gcf,'color','white');
% % 
mixtureModel = @(x,y) n1/(n1+n2+n3)*circ_bvmpdf(x,y,mu(1),nu(1),kappa1(1),kappa2(1),kappa3(1)) + ...
    n2/(n1+n2+n3)*circ_bvmpdf(x,y,mu(2),nu(2),kappa1(2),kappa2(2),kappa3(2)) + ...
    n3/(n1+n2+n3)*circ_bvmpdf(x,y,mu(2),nu(1),kappa1(2),kappa2(2),kappa3(2));
[xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
ppp = mixtureModel(xx,yy);
ppp = reshape(ppp,size(xx));
figure;contour(xx,yy,ppp,'ShowText','on');axis square;axis tight;
set(gcf,'color','white') 

[ params,posterior_probs, prior_probs] = mixture_of_bivariate_VM(alpha, 3);
est_mixtureModel = @(x,y) prior_probs(1)*circ_bvmpdf(x,y,params.mu(1),params.nu(1),params.kappa1(1),params.kappa2(1),params.kappa3(1)) + ...
    prior_probs(2)*circ_bvmpdf(x,y,params.mu(2),params.nu(2),params.kappa1(2),params.kappa2(2),params.kappa3(2)) + ...
    prior_probs(3)*circ_bvmpdf(x,y,params.mu(3),params.nu(3),params.kappa1(3),params.kappa2(3),params.kappa3(3));
[xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
ppp = est_mixtureModel(xx,yy);
ppp = reshape(ppp,size(xx));
figure;contour(xx,yy,ppp,'ShowText','on');axis square;axis tight;
set(gcf,'color','white') 

% alpha1 = circ_bvmrnd([mu(1) nu(1) kappa1(1) kappa2(1) 5], n1);
% fun = @(x,y,kp3) circ_bvmpdf(x,y,mu(1),nu(1),kappa1(1),kappa2(1),kp3);
% 
% [xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
% ppp = fun(xx,yy,5);
% ppp = reshape(ppp,size(xx));
% figure;contour(xx,yy,ppp,'ShowText','on');axis square;axis tight;
% set(gcf,'color','white') 
% figure;contour3(xx,yy,ppp)