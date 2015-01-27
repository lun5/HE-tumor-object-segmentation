

% $Author: ChrisMcCormick $    $Date: 2014/05/19 22:00:00 $    $Revision: 1.3 $
% Mixture of von Mises distribution with an uniform noise component

addpath(genpath(fullfile('Users','lun5','Research','github','HE-tumor-object-segmentation','toolboxes','CircStat2012a')));
%%======================================================
%% STEP 1a: Generate data from two 1D distributions.

mu1 = 0;      % Mean
kappa1 = 1;    % kappa
m1 = 100;      % Number of points

mu2 = pi/2;
kappa2 = 3;
m2 = 200;

% Generate the data.
X1 = circ_vmrnd(mu1,kappa1, m1);
X2 = circ_vmrnd(mu2,kappa2, m2);

% have to convert this into cartesian coordiates
X1_cart = [cos(X1) sin(X1)];
X2_cart = [cos(X2) sin(X2)];
X = [X1; X2];
X_cart = [X1_cart; X2_cart]; 
d = 2; % dimension
%%=====================================================
%% STEP 1b: Plot the data points and their pdfs.

x = -pi:0.1:pi;
y1 = circ_vmpdf(x, mu1, kappa1);
y2 = circ_vmpdf(x, mu2, kappa2);

figure;
plot(x, y1, 'b-');
hold on;
plot(x, y2, 'r-');
plot(X1, zeros(size(X1)), 'bx', 'markersize', 10);
plot(X2, zeros(size(X2)), 'rx', 'markersize', 10);
xlim([-pi pi]);

set(gcf,'color','white') % White background for the figure.

%% INPUT
% data
% number of components (excluding noise)
% 
%%====================================================
%% STEP 2: Choose initial values for the parameters.

% Set 'm' to the number of data points.
m = size(X, 1);

% Set 'k' to the number of clusters to find.
k = 2;

% this part here is where the mixture of von Mises happen

%% Initialization
% do k-means clustering to assign probabilites of component memberships to
% each of the n observations
% a-posteriori probabilities, results of spherical k-means
% don't have it so I will use linear k-means now
idx = kmeans(X, k);
posterior_probs = zeros(m,k);
for i = 1:k
    posterior_probs(idx == i,i) = 1;
end
likelihood = sum(posterior_probs,2);

%% Loop through M-step and E-step until convergence
maxiter = 1000;
% probability of each cluster -- prior
prior_probs = one(1,k)*(1/k);
mu_hat = zeros(1,k);
kappa_hat = zeros(1,k);
ep1 = 1e-4; % threshold for likelihood
ep2 = 1e-4; % threshold for parameters

for iter = 1: maxiter
    
    mu_hat_old = mu_hat;
    kappa_hat_old = kappa_hat;
    %% M-step
    for i = 1:k
        prior_probs(i) = mean(posterior_probs(:,i));
        unnormalized_mean = sum(repmat(posterior_probs(:,1),1,d).*X_cart);
        mu_hat(i) = unnormalized_mean/norm(unnormalized_mean);
        rho = norm(unnormalized_mean)/sum(posterior_probs(:,i));
        kappa_hat(i) = rho*(d - rho^2)/(1-rho^2);
    end
    
    %% E-step
    for i = 1:k
        posterior_probs(:,i) = prior_prob(i)*circ_vmpdf(X, mu_hat(i), kappa_hat(i));
    end
    
    likelihood_old = likelihood;
    likelihood = sum(posterior_probs,2);
    
    %% Stopping criteria
    llh_change = norm(abs(log(likelihood) - log(likelihood_old))) - ep1;
    mu_change = norm(abs(mu_hat - mu_hat_old)) - ep2;
    kappa_change = norm(abs(kappa_hat - kappa_hat_old)) - ep2;
    
    if llh_change || mu_change || kappa_change
        break;
    end
    
  
end

%%=====================================================
%% Plot the data points and their estimated pdfs.

x = -pi:0.1:pi;
y1 = circ_vmpdf(x, mu_hat(1), kappa_hat(1));
y2 = circ_vmpdf(x, mu_hat(2), kappa_hat(2));

% Plot over the existing figure, using black lines for the estimated pdfs.
plot(x, y1, 'k-');
plot(x, y2, 'k-');

hold off;

%%
% indeces = randperm(m);
% mu = zeros(1, k);
% for i = 1 : k
%     mu(i) = X(indeces(i));
% end
% 
% % Use the overal variance of the dataset as the initial variance for each cluster.
% kappa = ones(1, k) * sqrt(circ_var(X));
% 
% % Assign equal prior probabilities to each cluster.
% phi = ones(1, k) * (1 / (k+1));
% 
% %%===================================================
% %% STEP 3: Run Expectation Maximization
% maxiter = 1000;
% % Matrix to hold the probability that each data point belongs to each cluster.
% % One row per data point, one column per cluster.
% W = zeros(m, k+1);
% 
% % Loop until convergence.
% for iter = 1:maxiter
%     
%     fprintf('  EM Iteration %d\n', iter);
% 
%     %%===============================================
%     %% STEP 3a: Expectation
%     %
%     % Calculate the probability for each data point for each distribution.
%     
%     % Matrix to hold the pdf value for each every data point for every cluster.
%     % One row per data point, one column per cluster.
%     pdf = zeros(m, k+1);
%     
%     % For each cluster...
%     for j = 1 : k
%         
%         % Evaluate the Gaussian for all data points for cluster 'j'.
%         pdf(:, j) = circ_vmpdf(X, mu(j), kappa(j));
%     end
%     
%     pdf(:,k+1) = 1/(2*pi); % noise component
%     
%     % Multiply each pdf value by the prior probability for each cluster.
%     %    pdf  [m  x  k]
%     %    phi  [1  x  k]   
%     %  pdf_w  [m  x  k]
%     pdf_w = bsxfun(@times, pdf, phi);
%     
%     % Divide the weighted probabilities by the sum of weighted probabilities for each cluster.
%     %   sum(pdf_w, 2) -- sum over the clusters.
%     W = bsxfun(@rdivide, pdf_w, sum(pdf_w, 2));
%     
%     %%===============================================
%     %% STEP 3b: Maximization
%     %%
%     %% Calculate the probability for each data point for each distribution.
% 
%     % Store the previous means so we can check for convergence.
%     prevMu = mu;    
%     
%     % For each of the clusters...
%     for j = 1 : k
%     
%         % Calculate the prior probability for cluster 'j'.
%         phi(j) = mean(W(:, j));
%         
%         % Calculate the new mean for cluster 'j' by taking the weighted
%         % average of *all* data points.
%         % CHANGE THIS USING CIRCULAR 
%         mu(j) = circ_mean(X,W(:,j));
%     
%         % Calculate the variance for cluster 'j' by taking the weighted
%         % average of the squared differences from the mean for all data
%         % points.
%         % CHANGE THIS USING CIRCULAR 
%         variance = circ_var(X,W(:,j));
%         
%         % USE EQN 2 OF THE PAPER
%         % Calculate kappa by taking the square root of the variance.
%         kappa(j) = sqrt(variance);
%     end         
%     
%     phi(k+1) = 1 - sum(phi(1:k));
%     
%     % Check for convergence.
%     % Comparing floating point values for equality is generally a bad idea, but
%     % it seems to be working fine.
%     % CHANGE THE CONVERGENCE CRITERIA AS WELL
%     if (abs(mu-prevMu) < 1e-5) % should look at the log likelihood instead
%         break
%     end
% 
% % End of Expectation Maximization loop.    
% end
% 
% %%=====================================================
% %% STEP 4: Plot the data points and their estimated pdfs.
% 
% x = [0:0.1:30];
% y1 = gaussian1D(x, mu(1), kappa(1));
% y2 = gaussian1D(x, mu(2), kappa(2));
% 
% % Plot over the existing figure, using black lines for the estimated pdfs.
% plot(x, y1, 'k-');
% plot(x, y2, 'k-');
% 
