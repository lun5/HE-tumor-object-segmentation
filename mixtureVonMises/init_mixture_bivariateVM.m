% initialize the clusters
opts_kmeans = statset('Display','final');
numClusters = k;%6;
[idx,C] = kmeans(data,numClusters,'Distance','cityblock',...
    'Replicates',5,'Options',opts_kmeans);

figure;
colors = {'k','r','g','b','y','c'};
for cl = 1:numClusters
    plot(data(idx==cl,1),data(idx==cl,2),'.','MarkerFaceColor',colors{cl},'MarkerSize',12);
    hold on
end
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3)
%legend('Cluster 1','Cluster 2','Centroids',...
%       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off

for cl = 1: numClusters
    %Use page 171 on Statistics of Bivariate von Mises
    mu_hat(cl) = C(cl,1);
    nu_hat(cl) = C(cl,2);
    S1_bar = mean(sin(data(idx==cl,1) - mu_hat(cl)).^2);
    S2_bar = mean(sin(data(idx==cl,2) - nu_hat(cl)).^2);
    S12_bar = mean(sin(data(idx==cl,1) - mu_hat(cl)).*sin(data(idx==cl,2) - nu_hat(cl)));
    kappa1_hat(cl) = (S2_bar - S12_bar)./(S1_bar*S2_bar - S12_bar^2);
    kappa2_hat(cl) = (S1_bar - S12_bar)./(S1_bar*S2_bar - S12_bar^2);
    kappa3_hat(cl) = - S12_bar./(S1_bar*S2_bar - S12_bar^2);
end

%
%     numClusters = 3;
%     %alldata = [data(:,1); data(:,2)];
%     [ mu_hat_polar,~, kappa_hat,~, ~] = moVM([cos(data(:,1)) sin(data(:,1))],numClusters);
%     posterior_probs = zeros(numData,k + opts.noise);
%     init_params.theta_hat = mu_hat_polar;
%     init_params.kappa_hat = kappa_hat;
%     
%     % probability of each cluster -- prior
%     prior_probs = ones(1,k+opts.noise)*(1/(k+ opts.noise));
%     mu_hat = zeros(k,1);
%     nu_hat = zeros(k,1);
%     kappa1_hat = zeros(k,1);
%     kappa2_hat = zeros(k,1);    
%     kappa3_hat = zeros(k,1);  
%     % init_params is a struct with theta_hat and kappa_hat fields
%     [mean_sorted,ind] = sort(init_params.theta_hat,'descend');
%     kappa_sorted = init_params.kappa_hat(ind);
%     
%     if k > length(mean_sorted)*2
%         error('Something wrong with your inputs');
%     end
%     
%     % assign mu, nu, and all the kappa's
%     mu_hat(1:length(mean_sorted)) = mean_sorted(1);
%     mu_hat(length(mean_sorted)+1) = mean_sorted(2);
%     mu_hat(length(mean_sorted)+2) = mean_sorted(2);
%     mu_hat(k) = mean_sorted(3);
%     
%     nu_hat(1:length(mean_sorted)) = mean_sorted;
%     nu_hat(length(mean_sorted)+1) = mean_sorted(2);
%     nu_hat(length(mean_sorted)+2) = mean_sorted(3);
%     nu_hat(k) = mean_sorted(3);
%     
%     max_kappa = 100; min_kappa = 2; 
%     kappa1_hat(1:length(mean_sorted)) = min(max_kappa,kappa_sorted(1));
%     kappa1_hat(length(mean_sorted)+1) = min(max_kappa,kappa_sorted(2));
%     kappa1_hat(length(mean_sorted)+2) = min(max_kappa,kappa_sorted(2));
%     kappa1_hat(k) = kappa_sorted(3);
%     
%     kappa2_hat(1:length(mean_sorted)) = min(max_kappa,kappa_sorted);
%     kappa2_hat(length(mean_sorted)+1) = min(max_kappa,kappa_sorted(2));
%     kappa2_hat(length(mean_sorted)+2) = min(max_kappa,kappa_sorted(3));
%     kappa2_hat(k) = min(max_kappa,kappa_sorted(3));
%     
%     kappa3_choices = [-1 -0.5 0.5 1]; ind_ran = randi(4);
%     kappa3_hat(:) = kappa3_choices(ind_ran);%(kappa1_hat + kappa2_hat)/2