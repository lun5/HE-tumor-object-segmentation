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

function [ params,posterior_probs, prior_probs] = mixture_of_bivariate_VM(data, k, opts)
    
    opts_default.maxiter = 100;
    opts_default.eps = 1e-2; % threshold for likelihood convergence
    opts_default.noise = 0;
   
    if nargin <4
        opts = opts_default;
    elseif nargin <3
        error('Function needs at least 3 inputs: data, number of components, initial values for the parameters');
    end
    
    if ~exist('opts','var')
        opts = opts_default;
    end
    
    if ~isfield(opts,'maxiter')
        opts.maxiter = opts_default.maxiter;
    end
    
    if ~isfield(opts,'eps')
        opts.eps = opts_default.eps;
    end
    
    if ~isfield(opts,'noise');
        opts.noise = opts_default.noise;
    end
    
    % check that the data is in -pi and pi range
    if sum(abs(data(:)) > pi) > 0 || size(data,2) > 2
        error('Wrong input for the data: need to be angular variables')
    end

    % number of data points.
    numData = size(data, 1);
    
    %% STEP 1: Initialization
    numClusters = 3;%3;
    [ mu_hat_polar,~, kappa_hat,~, ~] = moVM([cos(data(:,1)) sin(data(:,1))],numClusters);
    posterior_probs = zeros(numData,k + opts.noise);
    init_params.theta_hat = mu_hat_polar;
    init_params.kappa_hat = kappa_hat;
    
    %% Loop through M-step and E-step until convergence
    % probability of each cluster -- prior
    prior_probs = ones(1,k+opts.noise)*(1/(k+ opts.noise));
    mu_hat = zeros(k,1);
    nu_hat = zeros(k,1);
    kappa1_hat = zeros(k,1);
    kappa2_hat = zeros(k,1);    
    kappa3_hat = zeros(k,1);  
    % init_params is a struct with theta_hat and kappa_hat fields
    [mean_sorted,ind] = sort(init_params.theta_hat,'descend');
    kappa_sorted = init_params.kappa_hat(ind);
    
    if k > length(mean_sorted)*2
        error('Something wrong with your inputs');
    end
    
    % assign mu, nu, and all the kappa's
    mu_hat(1:length(mean_sorted)) = mean_sorted(1);
    mu_hat(length(mean_sorted)+1) = mean_sorted(2);
    mu_hat(length(mean_sorted)+2) = mean_sorted(2);
    mu_hat(k) = mean_sorted(3);
    
    nu_hat(1:length(mean_sorted)) = mean_sorted;
    nu_hat(length(mean_sorted)+1) = mean_sorted(2);
    nu_hat(length(mean_sorted)+2) = mean_sorted(2);
    nu_hat(k) = mean_sorted(3);
    
    kappa1_hat(1:length(mean_sorted)) = kappa_sorted(1);
    kappa1_hat(length(mean_sorted)+1) = kappa_sorted(2);
    kappa1_hat(length(mean_sorted)+2) = kappa_sorted(2);
    kappa1_hat(k) = kappa_sorted(3);
    
    kappa2_hat(1:length(mean_sorted)) = kappa_sorted;
    kappa2_hat(length(mean_sorted)+1) = kappa_sorted(2);
    kappa2_hat(length(mean_sorted)+2) = kappa_sorted(2);
    kappa2_hat(k) = kappa_sorted(3);
    
    kappa3_hat(:) = 1;%(kappa1_hat + kappa2_hat)/2;
    weighted_logllh = zeros(k,1);
    
    fmincon_opts = optimset('PlotFcns',@optimplotfval,'Display','iter',...
        'MaxIter',500,'TolFun',1e-5,'TolX',1e-5);
    % linear constraints for fmincon
    %A = [0 0 -1 0 1; 0 0 0 -1 1]; b = [0 0];
    % non linear constraint defined in the end
    % upper and lower bounds
    offset_mean = 0.5; offset_conc = 3;
    lb = [mu_hat - offset_mean nu_hat - offset_mean max(kappa1_hat -offset_conc,0) max(kappa2_hat -offset_conc,0) ones(size(mu_hat))]; 
    ub = [mu_hat + offset_mean nu_hat + offset_mean ones(size(mu_hat))*max(50,max(kappa1_hat)) ones(size(mu_hat))*max(50,max(kappa2_hat)) ones(size(mu_hat))*2];
%     lb = [mu_hat - offset_mean nu_hat - offset_mean ones(size(mu_hat))*10 ones(size(mu_hat))*10 ones(size(mu_hat))]; 
%     ub = [mu_hat + offset_mean nu_hat + offset_mean ones(size(mu_hat))*50 ones(size(mu_hat))*50 ones(size(mu_hat))*2]; 

for iter = 1: opts.maxiter
    mu_hat_old = mu_hat;
    nu_hat_old = nu_hat;
    kappa1_hat_old = kappa1_hat;
    kappa2_hat_old = kappa2_hat;
    kappa3_hat_old = kappa3_hat;
    weighted_logllh_old = weighted_logllh;
    %% STEP 1: E-step
    for i = 1:k
        posterior_probs(:,i) = prior_probs(i)*circ_bvmpdf(data(:,1),data(:,2),...
            mu_hat(i), nu_hat(i), kappa1_hat(i),  kappa2_hat(i), kappa3_hat(i));
    end
    if opts.noise
        posterior_probs(:,k+1) = prior_probs(k+1)*repmat(1/(2*pi)^2,numData,1);
    end
    
    posterior_probs = posterior_probs./repmat(sum(posterior_probs,2),1,k + opts.noise);
    %% STEP 2: M-step
    for i = 1:k
        prior_probs(i) = mean(posterior_probs(:,i));
        
%         [param_best, funcval_final, exitflag] = fminsearch(@(x) ...
%             - Loglikelihood(posterior_probs(:,i),data(:,1),data(:,2),[x kappa3_hat(i)]),...
%             [mu_hat(i),nu_hat(i),kappa1_hat(i),kappa2_hat(i)], fminsearch_opts);
        [param_best, funcval_final, exitflag] = fmincon(@(x) ...
            - Loglikelihood(posterior_probs(:,i),data(:,1),data(:,2),x),...
            [mu_hat(i),nu_hat(i),kappa1_hat(i),kappa2_hat(i),kappa3_hat(i)],...
            [],[],[],[],lb(i,:),ub(i,:),@confun,fmincon_opts);

        if exitflag ~=1 % run it twice at most
            [param_best, funcval_final, exitflag] = fmincon(@(x) ...
            - Loglikelihood(posterior_probs(:,i),data(:,1),data(:,2),x),...
            param_best,[],[],[],[],lb(i,:),ub(i,:),@confun,fmincon_opts);
%              [param_best, funcval_final, exitflag] = fminsearch(@(x) ...
%             - Loglikelihood(posterior_probs(:,i),data(:,1),data(:,2),[x kappa3_hat(i)]),...
%             param_best, fminsearch_opts);
        end        
        weighted_logllh(i) = funcval_final;
        
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
    llh_change = norm(abs(weighted_logllh - weighted_logllh_old));
    mu_change = norm(abs(mu_hat - mu_hat_old));
    nu_change = norm(abs(nu_hat - nu_hat_old));
    kappa1_change = norm(abs(kappa1_hat - kappa1_hat_old));
    kappa2_change = norm(abs(kappa2_hat - kappa2_hat_old));
    kappa3_change = norm(abs(kappa3_hat - kappa3_hat_old));
    
    if llh_change  < opts.eps|| mu_change  < opts.eps || nu_change  < opts.eps ...
            || kappa1_change < opts.eps || kappa2_change  < opts.eps %|| kappa3_change  < opts.eps
        break;
    end
  
end

params.mu = mu_hat;
params.nu = nu_hat;
params.kappa1 = kappa1_hat;
params.kappa2 = kappa2_hat;
params.kappa3 = kappa3_hat;

if iter == opts.maxiter
    sprintf('The algorithm does not converge at maxiter %d',opts.maxiter)
end
  
end

function [LLH] = Loglikelihood(pij,phi, psi,params)
    mu = params(1); nu = params(2); kappa1 = params(3); kappa2 = params(4); kappa3 = params(5);
    H = kappa1*cos(phi-mu) + kappa2*cos(psi - nu) - kappa3*cos(phi-mu-psi+nu);
    fun = @(x, nu, kappa1, kappa2, kappa3) 2*pi*besseli(0,sqrt(kappa1.^2+kappa3.^2 ...
    -2*kappa1.*kappa3.*cos(x - nu))).*exp(kappa2.*cos(x-nu));
    Cc_inv = integral((@(x)fun(x, nu, kappa1, kappa2, kappa3)),0,2*pi);
    ind = pij > 1e-5;
    LLH = length(phi(ind))*log(Cc_inv^-1) + sum(H(ind)) + sum(log(pij(ind))); % log likelihood
end

function [c, ceq] = confun(params)
kappa1 = params(3); kappa2 = params(4); kappa3 = params(5);
c = [kappa3 - kappa1*kappa2/(kappa1+kappa2)];%;kappa3 - kappa1;kappa3 - kappa2];
ceq = [];
end
