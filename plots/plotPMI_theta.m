% INPUT: sample, the parameters of the mixture model
%function plotPMI_theta
    Fsym = [F; F(:,2) F(:,1)];
    %figure; ndhist(Fsym(:,1),Fsym(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
    figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
    xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
    set(gcf,'color','white') 
    % save this into a file
    
    % mixture model 
%     est_mixtureModel = @(x,y) prior_probs(1)*circ_bvmpdf(x,y,params.mu(1),params.nu(1),params.kappa1(1),params.kappa2(1),params.kappa3(1)) + ...
%     prior_probs(2)*circ_bvmpdf(x,y,params.mu(2),params.nu(2),params.kappa1(2),params.kappa2(2),params.kappa3(2)) + ...
%     prior_probs(3)*circ_bvmpdf(x,y,params.mu(3),params.nu(3),params.kappa1(3),params.kappa2(3),params.kappa3(3)) + ...
%     prior_probs(4)*circ_bvmpdf(x,y,params.mu(4),params.nu(4),params.kappa1(4),params.kappa2(4),params.kappa3(4)) + ...
%     prior_probs(5)*circ_bvmpdf(x,y,params.mu(5),params.nu(5),params.kappa1(5),params.kappa2(5),params.kappa3(5)) + ...
%     prior_probs(6)*circ_bvmpdf(x,y,params.mu(6),params.nu(6),params.kappa1(6),params.kappa2(6),params.kappa3(6));

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
    yy = - yy; 
%     for i = 1:6
%         ppp = prior_probs(i) * circ_bvmpdf(xx,yy,params.mu(i),params.nu(i),...
%             params.kappa1(i),params.kappa2(i),params.kappa3(i));
%         ppp = reshape(ppp,size(xx)); %ppp = (ppp+ppp')/2;
%         figure;contourf(xx,yy,ppp);
%         axis square;axis tight;
%         set(gcf,'color','white');colorbar;
%         xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
%     end
%     
    ppp = est_mixtureModel(xx,yy);
    ppp = reshape(ppp,size(xx)); %ppp = (ppp+ppp')/2;
    numContours = 10;
    figure;mesh(xx,yy,ppp);%
    %figure; contourf(xx,yy,ppp,numContours);colorbar;
    axis square;axis tight;
    set(gcf,'color','white');
    xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
       
%     figure;contourf(xx,yy,log(ppp),numContours);
%     axis square;axis tight;
%     set(gcf,'color','white');colorbar;
%     xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
    
    [pmi,pJoint,pProd] = evalPMI_theta([xx(:),yy(:)], mixture_params,opts);
    ppp = reshape(pmi,size(xx)); %ppp = (ppp + ppp')./2;
    figure;contourf(xx,yy,ppp,numContours);axis square;axis tight;
    set(gcf,'color','white');colorbar;
    xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
    %caxis([0 0.5]); 
%     figure;contourf(xx,yy,log(ppp),numContours);axis square;axis tight;
%     set(gcf,'color','white');colorbar;
%     xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
    %xlabel('\phi'); ylabel('\psi');set(gca,'FontSize',16);    
%end
