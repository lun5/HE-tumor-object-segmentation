% %% Interactive selection of red, green, blue circles on the zebra
% figure; imshow(I./255); hold on;
% 
% [rowSub,colSub] = ginput;
% rowSub = round(rowSub); colSub = round(colSub); 
% c_vecs = {'r','g','b'};
% % shape inserter for 3 combinations: 
% % pink-pink: red circle, purple purple: white circle, pink-purple: green
% % circle
% 
% for i =1:floor(length(colSub)/2)
%     coord1 = [rowSub(2*(i-1)+1),colSub(2*(i-1)+1)];
%     coord2 = [rowSub(2*i),colSub(2*i)];
%     r_shape = 30; % radius of the shape drawn
%     plot([coord1(1),coord2(1)],[coord1(2),coord2(2)],'o','MarkerSize',5,...
%         'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});    
%         rectangle('Position',[(coord1(1) + coord2(1))/2 - r_shape ...
%             (coord1(2) + coord2(2))/2 - r_shape 2*r_shape 2*r_shape],...,
%             'LineWidth', 5, 'EdgeColor',c_vecs{i})
% end
% 
% hold off;
% 
% f_maps_cur = f_maps{1};
% linearInd = sub2ind(size(f_maps_cur), colSub, rowSub);
% feats = f_maps_cur(linearInd);

%% INPUT: sample, the parameters of the mixture model
%function plotPMI_theta
Fsym = [F; F(:,2) F(:,1)];
%figure; ndhist(Fsym(:,1),Fsym(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 
 
% mixture model 
if opts.model_half_space_only
    est_mixtureModel = @(x,y) prior_probs(1)*circ_bvmpdf(x,y,params.mu(1),params.nu(1),params.kappa1(1),params.kappa2(1),params.kappa3(1)) + ...
    prior_probs(2)*circ_bvmpdf(x,y,params.mu(2),params.nu(2),params.kappa1(2),params.kappa2(2),params.kappa3(2)) + ...
    prior_probs(3)*circ_bvmpdf(x,y,params.mu(3),params.nu(3),params.kappa1(3),params.kappa2(3),params.kappa3(3)) + ...
    prior_probs(4)*circ_bvmpdf(x,y,params.mu(4),params.nu(4),params.kappa1(4),params.kappa2(4),params.kappa3(4)) + ...
    prior_probs(5)*circ_bvmpdf(x,y,params.mu(5),params.nu(5),params.kappa1(5),params.kappa2(5),params.kappa3(5)) + ...
    prior_probs(6)*circ_bvmpdf(x,y,params.mu(6),params.nu(6),params.kappa1(6),params.kappa2(6),params.kappa3(6));
else
    est_mixtureModel = @(x,y) prior_probs(1)*circ_bvmpdf(x,y,params.mu(1),params.nu(1),params.kappa1(1),params.kappa2(1),params.kappa3(1)) + ...
     prior_probs(2)*circ_bvmpdf(x,y,params.mu(2),params.nu(2),params.kappa1(2),params.kappa2(2),params.kappa3(2)) + ...
     prior_probs(3)*circ_bvmpdf(x,y,params.mu(3),params.nu(3),params.kappa1(3),params.kappa2(3),params.kappa3(3)) + ...
     prior_probs(4)*circ_bvmpdf(x,y,params.mu(4),params.nu(4),params.kappa1(4),params.kappa2(4),params.kappa3(4)) + ...
     prior_probs(5)*circ_bvmpdf(x,y,params.mu(5),params.nu(5),params.kappa1(5),params.kappa2(5),params.kappa3(5)) + ...
     prior_probs(6)*circ_bvmpdf(x,y,params.mu(6),params.nu(6),params.kappa1(6),params.kappa2(6),params.kappa3(6)) + ...
     prior_probs(7)*circ_bvmpdf(x,y,params.mu(7),params.nu(7),params.kappa1(7),params.kappa2(7),params.kappa3(7)) + ...
     prior_probs(8)*circ_bvmpdf(x,y,params.mu(8),params.nu(8),params.kappa1(8),params.kappa2(8),params.kappa3(8)) + ...
     prior_probs(9)*circ_bvmpdf(x,y,params.mu(9),params.nu(9),params.kappa1(9),params.kappa2(9),params.kappa3(9));
end
%%
[xx,yy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
%     for i = 1:6
%         ppp = prior_probs(i) * circ_bvmpdf(xx,yy,params.mu(i),params.nu(i),...
%             params.kappa1(i),params.kappa2(i),params.kappa3(i));
%         ppp = reshape(ppp,size(xx)); %ppp = (ppp+ppp')/2;
%         figure;contourf(xx,yy,ppp);
%         axis square;axis tight;
%         set(gcf,'color','white');colorbar;
%         xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
%     end
     
ppp = est_mixtureModel(xx,yy);
ppp = reshape(ppp,size(xx)); 
if opts.model_half_space_only; ppp = (ppp+ppp')/2; end
numContours = 10;
figure;mesh(xx,yy,ppp);%
figure; contourf(xx,yy,ppp,numContours);colorbar;%caxis([0 0.4])
axis square;axis tight;
set(gcf,'color','white');
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
% hold on;
% for i =1:floor(length(feats)/2)
%     coord = [feats(2*(i-1)+1) feats(2*i)]; % coordinate of points picked interactively
%     r_shape = 0.3;
%     plot(coord(1),coord(2),'o','MarkerSize',4,...
%         'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
%         rectangle('Position', [ coord(1) - r_shape coord(2) - r_shape ...
%             2*r_shape 2*r_shape], 'LineWidth', 2, 'EdgeColor',c_vecs{i});
% end
% hold off;   

%%
[pmi,pJoint,pProd] = evalPMI_theta([xx(:),yy(:)], mixture_params{1},opts);
%[pmi,pJoint,pProd] = approxPMI_theta([xx(:),yy(:)], mixture_params,opts);

ppp = reshape(pmi,size(xx));%ppp = reshape(pJoint,size(xx));ppp = reshape(pProd,size(xx));
if opts.model_half_space_only; ppp = (ppp+ppp')/2; end
figure;contourf(xx,yy,ppp,numContours);axis square;axis tight;
set(gcf,'color','white');colorbar;
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',20);
    
% hold on;
% for i =1:floor(length(feats)/2)
%     coord = [feats(2*(i-1)+1) feats(2*i)]; % coordinate of points picked interactively
%     r_shape = 0.3;
%     plot(coord(1),coord(2),'o','MarkerSize',4,...
%         'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
%         rectangle('Position', [ coord(1) - r_shape coord(2) - r_shape ...
%             2*r_shape 2*r_shape], 'LineWidth', 2, 'EdgeColor',c_vecs{i});
% end
% hold off;

% pJoint
ppp = reshape(log(pJoint),size(xx));
if opts.model_half_space_only; ppp = (ppp+ppp')/2; end
figure;contourf(xx,yy,ppp,numContours);axis square;axis tight;
set(gcf,'color','white');colorbar;
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',20);
%     
% hold on;
% for i =1:floor(length(feats)/2)
%     coord = [feats(2*(i-1)+1) feats(2*i)]; % coordinate of points picked interactively
%     r_shape = 0.3;
%     plot(coord(1),coord(2),'o','MarkerSize',4,...
%         'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
%         rectangle('Position', [ coord(1) - r_shape coord(2) - r_shape ...
%             2*r_shape 2*r_shape], 'LineWidth', 4, 'EdgeColor',c_vecs{i});
% end
% hold off;
% 
% pMarge
ppp = reshape(pProd,size(xx));
if opts.model_half_space_only; ppp = (ppp+ppp')/2; end
figure;contourf(xx,yy,ppp,numContours);axis square;axis tight;
set(gcf,'color','white');colorbar;
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
    
% hold on;
% for i =1:floor(length(feats)/2)
%     coord = [feats(2*(i-1)+1) feats(2*i)]; % coordinate of points picked interactively
%     r_shape = 0.3;
%     plot(coord(1),coord(2),'o','MarkerSize',4,...
%         'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
%         rectangle('Position', [ coord(1) - r_shape coord(2) - r_shape ...
%             2*r_shape 2*r_shape], 'LineWidth', 2, 'EdgeColor',c_vecs{i});
% end
% hold off;
% 



