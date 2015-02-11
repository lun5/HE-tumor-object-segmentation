%function [ params,posterior_probs, prior_probs] = moVM2D(data,k,init_params, opts)
% 
%  data          - n data points x 2 --> bivariate von Mises:
%                - number of data points. Range -pi to pi
%                    
%  k             - number of clusters to find.
%  init_param    - initial parameters from 
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

function [ params,posterior_probs, prior_probs] = moVM2D(data,k,init_params, opts)
    
    opts_default.maxiter = 100;
    opts_default.eps1 = 1e-2; % threshold for likelihood convergence
    opts_default.eps2 = 1e-2; % threshold for parameter convergence
    opts_default.noise = 1;
    
    if nargin <4
        opts = opts_default;
    elseif nargin <3
        error('Function needs at least 3 inputs: data, number of components, initial values for the parameters');
    end
    
    if ~exist('opts','var')
        opts = opts_default;
    end
    
    if ~exist('opts.maxiter','var')
        opts.maxiter = opts_default.maxiter;
    end
    
    if ~exist('opts.eps1','var');
        opts.eps1 = opts_default.eps1;
    end

    if ~exist('opts.eps2','var');
        opts.eps2 = opts_default.eps2;
    end
    
    if ~exist('opts.noise','var');
        opts.noise = opts_default.noise;
    end
    
    % check that the data is in -pi and pi range
    if abs(data(:)) > pi || size(data,2) > 2
        error('Wrong input for the data: need to be angular variables')
    end

    % Set 'm' to the number of data points.
    numData = size(data, 1);
    
    %% STEP 1: Initialization
    posterior_probs = zeros(numData,k + opts.noise);
    
    %% Loop through M-step and E-step until convergence
    % probability of each cluster -- prior
    prior_probs = ones(1,k+opts.noise)*(1/(k+ opts.noise));
    mu_hat = zeros(k,1);
    nu_hat = zeros(k,1);
    kappa1_hat = zeros(k,1);
    kappa2_hat = zeros(k,1);
    kappa3_hat = zeros(k,1);
    
    % init_params is a struct with theta_hat and kappa_hat fields
    [mean_sorted,ind] = sort(init_params.theta_hat);
    kappa_sorted = init_params.kappa_hat(ind);
    
    if k > length(mean_sorted)*2
        error('Something wrong with your inputs');
    end
    
    % assign mu, nu, and all the kappa's
    mu_hat(1:length(mean_sorted)) = mean_sorted(1);
    mu_hat(length(mean_sorted)+1:length(mean_sorted)+2) = mean_sorted(2);
    mu_hat(k) = mean_sorted(3);
    
    nu_hat(1:length(mean_sorted)) = mean_sorted;
    nu_hat(length(mean_sorted)+1:length(mean_sorted)+2) = mean_sorted(2:3);
    nu_hat(k) = mean_sorted(3);
    
    kappa1_hat(1:length(mean_sorted)) = kappa_sorted(1);
    kappa1_hat(length(mean_sorted)+1:length(mean_sorted)+2) = kappa_sorted(2);
    kappa1_hat(k) = kappa_sorted(3);
    
    kappa2_hat(1:length(mean_sorted)) = kappa2_sorted;
    kappa2_hat(length(mean_sorted)+1:length(mean_sorted)+2) = kappa2_sorted(2:3);
    kappa2_hat(k) = kappa2_sorted(3);
    
    kappa3_hat = (kappa1_hat + kappa2_hat)/2;
    
for iter = 1: opts.maxiter
    
    mu_hat_old = mu_hat;
    nu_hat_old = nu_hat;
    kappa1_hat_old = kappa1_hat;
    kappa2_hat_old = kappa2_hat;
    kappa3_hat_old = kappa3_hat;
    
    %% STEP 1: E-step
    posterior_probs_old = posterior_probs;
    for i = 1:k
        posterior_probs(:,i) = prior_probs(i)*circ_bvmpdf(data(:,1),data(:,2),...
            mu_hat(i), nu_hat(i), kappa1_hat(i),  kappa2_hat(i), kappa3_hat(i));
    end
    if opts.noise
        posterior_probs(:,k+1) = prior_probs(k+1)*repmat(1/(2*pi)^2,numData,1);
    end
    
    posterior_probs = posterior_probs./repmat(sum(posterior_probs,2),1,k + opts.noise);
    fminsearch_opts.Display = 'off';%'iter';
    fminsearch_opts.MaxIter = 20;
    %% STEP 2: M-step
    for i = 1:k
        prior_probs(i) = mean(posterior_probs(:,i));
        f = @(mu,nu,kappa1,kappa2,kappa3) prod(LLikelihood(posterior_probs(:,i), data(:,1),...
            data(:,2),mu, nu, kappa1, kappa2, kappa3));
        [params] = fminsearch(@(bw) f(bw,p,F_val), bw(1:size(bw,1)/2), fminsearch_opts);
        [param_best, funcval_final, exitflag] = fminsearch(@(x) f(mu,nu,kappa1,kappa2,kappa3),...
            [mu_hat(i),nu_hat(i),kappa1_hat(i),kappa2_hat(i),kappa3_hat(i)], fminsearch_opts);
        %% save the funcval somewhere!!!
        mu_hat(i) = param_best(1);
        nu_hat(i) = param_best(2);
        kappa1_hat(i) = param_best(3);
        kappa2_hat(i) = param_best(4);
        kappa3_hat(i) = param_best(5);
    end
    
    % rescale the uniform noise if it goes above 10%
    if sum(prior_probs(1:k)) < 0.8
        prior_probs(1:k) = prior_probs(1:k)*(0.8+0.1*rand)/sum(prior_probs(1:k));
    end
    
    if opts.noise
        prior_probs(k+1) = 1 - sum(prior_probs(1:k));
    end
        
    %% Stopping criteria
    % There is one very concentrated cluster of white. If the concentration
    % of this cluster is greater than 600 then we will stop the algorithm
%     if max(kappa_hat) > 300
%         break;
%         % or maybe use spkmeans at this point
%     end

    % else look at the change in posterior and parameters
    %llh_change = norm(abs(log(posterior_probs+1e-10) - log(posterior_probs_old+1e-10)));
    llh_change = norm(abs(posterior_probs - posterior_probs_old));
    mu_change = norm(abs(mu_hat_polar - mu_hat_old));
    kappa_change = norm(abs(kappa_hat - kappa_hat_old));
    
    if llh_change  < opts.eps1|| (mu_change  < opts.eps2 || kappa_change  < opts.eps2)
        break;
    end
  
end

if iter == opts.maxiter
    sprintf('The algorithm does not converge at maxiter %d',opts.maxiter)
end
  
end

function [H] = LLikelihood(pij, phi, psi,mu, nu, kappa1, kappa2, kappa3)
    H = pij.*circ_bvmpdf(phi,psi,mu,nu,kappa1,kappa2,kappa3);
end