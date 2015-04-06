% testing the relevant of kappa3
n = 500;
mu = 2; nu = 1;
kappa1 = 5; kappa2 = 2;
kappa3_vec = [-0.5 0 0.5];
x = linspace(1e-6-pi,pi); y = linspace(1e-6-pi,pi);
[xx,yy] = meshgrid(x,y);
params = zeros(1,5);
for i = 1:length(kappa3_vec)
    kappa3 = kappa3_vec(i);
    params(1) = mu; params(2) = nu;
    params(3) = kappa1; params(4) = kappa2; params(5) = kappa3;
%     [phi,psi] = circ_bvmrnd(params, n);
%     %histogram
%     figure; ndhist(phi,psi,'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
%     xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
%     set(gcf,'color','white')
    % contour plot
    zz = circ_bvmpdf(xx,yy,mu,nu,kappa1,kappa2,kappa3);
    zz = reshape(zz,size(xx));
    figure; contour(xx,yy,zz,10)
    xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
    set(gcf,'color','white')
end
    