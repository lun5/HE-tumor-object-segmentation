%addpath(genpath(pwd));
mu = [0 2];
nu = [0 2];
kappa1 = [1 1];
kappa2 = [1 1];
kappa3 = [-0.5 -0.5];
n1 = 100; n2 = 100;
% 
% alpha1 = circ_bvmrnd([mu(1) nu(1) kappa1(1) kappa2(1) kappa3(1)], n1);
% alpha2 = circ_bvmrnd([mu(2) nu(2) kappa1(2) kappa2(2) kappa3(2)], n2);
% alpha = [alpha1; alpha2];
% figure; ndhist(alpha(:,1),alpha(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
% xlim([-pi pi]); ylim([-pi pi]); xlabel('phi');ylabel('psi');set(gcf,'color','white') 
% 
% mixtureModel = @(x,y) n1/(n1+n2)*circ_bvmpdf(x,y,mu(1),nu(1),kappa1(1),kappa2(1),kappa3(1)) + ...
%     n1/(n1+n2)*circ_bvmpdf(x,y,mu(2),nu(2),kappa1(2),kappa2(2),kappa3(2));
% [xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
% ppp = mixtureModel(xx,yy);
% ppp = reshape(ppp,size(xx));
% figure;contour(xx,yy,ppp,'ShowText','on');axis square;axis tight;
% set(gcf,'color','white') 

alpha1 = circ_bvmrnd([mu(1) nu(1) kappa1(1) kappa2(1) 5], n1);
fun = @(x,y,kp3) circ_bvmpdf(x,y,mu(1),nu(1),kappa1(1),kappa2(1),kp3);

[xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
ppp = fun(xx,yy,5);
ppp = reshape(ppp,size(xx));
figure;contour(xx,yy,ppp,'ShowText','on');axis square;axis tight;
set(gcf,'color','white') 
figure;contour3(xx,yy,ppp)
