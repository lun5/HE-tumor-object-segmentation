%% function [ mu_hat_polar, mu_hat_cart, kappa_hat,posterior_probs] = moVM(X_cart,k,opts)
% 
% INPUTS ONLY IN 2D!!!!!
%  X             - data in cartesian space of size mxd. d: dimension, m:
%                - number of data points. Will convert to polar
%                    
%  k             - number of clusters to find.
%  opts          - parameter settings (maxiter, initialization)
%
% OUTPUTS
%  mu_hat_polar  - estimated means in polar space
%  kappa_hat     - estimated concentration parameters
%  mu_hat_cart   - estimated means in cartesian space
%  posterior_probs - component membership P(in kth component| data i)
%
% -------------------------------------------------------------------------
% HE segmentation toolbox
% Luong Nguyen, 2014 [lun5@pitt.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

function [ mu_hat_polar, mu_hat_cart, kappa_hat,posterior_probs, prior_probs] = moVM(X_cart,k,opts)
    
    opts_default.maxiter = 1000;
    opts_default.init = 'lkmeans'; % initialize by linear kmeans clustering
    opts_default.eps1 = 1e-4; % threshold for likelihood convergence
    opts_default.eps2 = 1e-5; % threshold for parameter convergence
    
    if nargin <3
        opts = opts_default;
    elseif nargin <2
        error('Function needs at least 2 inputs: data, number of components');
    end
    
    if ~exist('opts.maxiter','var')
        opts.maxiter = opts_default.maxiter;
    end
    
    if ~exist('opts.init','var');
        opts.init = opts_default.init;
    end
    
    if ~exist('opts.eps1','var');
        opts.eps1 = opts_default.eps1;
    end

    if ~exist('opts.eps2','var');
        opts.eps1 = opts_default.eps2;
    end

    % Set 'm' to the number of data points.
    m = size(X_cart, 1);
    % set 'd' to the dimension
    d = size(X_cart,2); % for now it's just 2D
    X_polar = atan2(X_cart(:,2),X_cart(:,1));

    %% STEP 1: Initialization
    % do k-means clustering to assign probabilites of component memberships to
    % each of the n observations
    % a-posteriori probabilities, results of spherical k-means
    % don't have it so I will use linear k-means now
    
%     posterior_probs = zeros(m,k+1);
%     if strcmp(opts.init,'lkmeans')
%         idx = kmeans(X_polar, k+1);
%         for i = 1:k
%            posterior_probs(idx == i,i) = 1;
%         end
%         posterior_probs(:,k+1) = 0; % uniform noise
%     else
%        posterior_probs(:) = 1/(k+1); 
%     end
    posterior_probs = zeros(m,k);
    if strcmp(opts.init,'lkmeans')
        idx = kmeans(X_polar, k);
        for i = 1:k
           posterior_probs(idx == i,i) = 1;
        end
    else
       posterior_probs(:) = 1/(k); 
    end

    %% Loop through M-step and E-step until convergence
    % probability of each cluster -- prior
    %prior_probs = ones(1,k+1)*(1/(k+1));
    prior_probs = ones(1,k)*(1/(k));    
    mu_hat_cart = zeros(d,k);
    mu_hat_polar = zeros(1,k);
    kappa_hat = zeros(1,k);
    
for iter = 1: opts.maxiter
    
    mu_hat_old = mu_hat_polar;
    kappa_hat_old = kappa_hat;
    %% STEP 2: M-step
    for i = 1:k
        prior_probs(i) = mean(posterior_probs(:,i));
        unnormalized_mean = sum(repmat(posterior_probs(:,i),1,d).*X_cart);
        mu_hat_cart(:,i) = unnormalized_mean'/norm(unnormalized_mean);
        mu_hat_polar(i) = atan2(mu_hat_cart(2,i),mu_hat_cart(1,i));
        rho = norm(unnormalized_mean)/sum(posterior_probs(:,i));
        kappa_hat(i) = rho*(d - rho^2)/(1-rho^2);
    end
    %prior_probs(k+1) = 1 - sum(prior_probs(1:k));
    
    %% STEP 1: E-step
    posterior_probs_old = posterior_probs;
    for i = 1:k
        posterior_probs(:,i) = prior_probs(i)*circ_vmpdf(X_polar, mu_hat_polar(i), kappa_hat(i));
    end
    %posterior_probs(:,k+1) = prior_probs(k+1)*repmat(1/(2*pi),m,1);
    
    % normalize posterior_probs such that sum of prior probs = 1
    %posterior_probs = posterior_probs./repmat(sum(posterior_probs,2),1,k+1);
    posterior_probs = posterior_probs./repmat(sum(posterior_probs,2),1,k);
    
    %% Stopping criteria
    %llh_change = norm(abs(log(posterior_probs+1e-10) - log(posterior_probs_old+1e-10)));
    llh_change = norm(abs(posterior_probs - posterior_probs_old));
    mu_change = norm(abs(mu_hat_polar - mu_hat_old));
    kappa_change = norm(abs(kappa_hat - kappa_hat_old));
    
    if llh_change  < opts.eps1 && (mu_change  < opts.eps2 && kappa_change  < opts.eps2)
        break;
    end
  
end

if iter == opts.maxiter
    disp('The algorithm does not converge at maxiter')
end
  
end