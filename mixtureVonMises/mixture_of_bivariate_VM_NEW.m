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

function [ params,posterior_probs, prior_probs] = mixture_of_bivariate_VM_NEW(data, k, opts)
    
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
    mu_hat = zeros(k,1); nu_hat = zeros(k,1); 
    kappa1_hat = zeros(k,1); kappa2_hat =  zeros(k,1); kappa3_hat = zeros(k,1);
    init_mixture_bivariateVM;  
    mult = 1.005; %for the annealing step
    threshold_kappa = mean(kappa1_hat);
%     numClusters = k;
%     opts_kmeans = statset('Display','final');
%     [idx,C] = kmeans(data,numClusters,'Distance','cityblock',...
%     'Replicates',5,'Options',opts_kmeans);
%     mu_hat = zeros(k,1);nu_hat = zeros(k,1);
%     kappa1_hat = zeros(k,1);kappa2_hat = zeros(k,1);kappa3_hat = zeros(k,1);
%     for cl = 1: numClusters
%         %Use page 171 on Statistics of Bivariate von Mises
%         mu_hat(cl) = C(cl,1);
%         nu_hat(cl) = C(cl,2);
%         S1_bar = mean(sin(data(idx==cl,1) - mu_hat(cl)).^2);
%         S2_bar = mean(sin(data(idx==cl,2) - nu_hat(cl)).^2);
%         S12_bar = mean(sin(data(idx==cl,1) - mu_hat(cl)).*sin(data(idx==cl,2) - nu_hat(cl)));
%         kappa1_hat(cl) = S2_bar./(S1_bar*S2_bar - S12_bar^2);
%         kappa2_hat(cl) = S1_bar./(S1_bar*S2_bar - S12_bar^2);
%         kappa3_hat(cl) = S12_bar./(S1_bar*S2_bar - S12_bar^2);
%     end
    
    weighted_logllh = zeros(k,1);
    posterior_probs = zeros(numData,k + opts.noise);
    % probability of each cluster -- prior
    prior_probs = ones(1,k+opts.noise)*(1/(k+ opts.noise));
    fmincon_opts = optimset('PlotFcns',@optimplotfval,'Display','iter',...
        'MaxIter',500);%,'TolFun',1e-5,'TolX',1e-5);
    %fmincon_opts = optimset('Display','on',...
    %     'MaxIter',500,'Algorithm','sqp');%,'TolFun',1e-5,'TolX',1e-5);
    % linear constraints for fmincon
    % non linear constraint defined in the end
    offset_mean = 0.2; 
    max_kappa = 100; min_kappa = 1;
    ub_kappa12 = max_kappa;% max(max_kappa,min(max_kappa,max(kappa_sorted))+ offset_conc);
    %lb = [mu_hat - offset_mean nu_hat - offset_mean min(kappa1_hat,min_kappa) min(kappa2_hat,min_kappa) ones(size(mu_hat))*(-2)];
    %ub = [mu_hat + offset_mean nu_hat + offset_mean ones(size(mu_hat))*ub_kappa12 ones(size(mu_hat))*ub_kappa12 ones(size(mu_hat))*2];
    lb = [mu_hat - offset_mean nu_hat - offset_mean min(kappa1_hat,min_kappa) min(kappa2_hat,min_kappa)]; % ones(size(mu_hat))*(0)];
    ub = [mu_hat + offset_mean nu_hat + offset_mean ones(size(mu_hat))*ub_kappa12 ones(size(mu_hat))*ub_kappa12]; % ones(size(mu_hat))*0]; 
    % should we use the observation that the shape is most elipse when
    
for iter = 1: opts.maxiter
    mu_hat_old = mu_hat;
    nu_hat_old = nu_hat;
    kappa1_hat_old = kappa1_hat;
    kappa2_hat_old = kappa2_hat;
    kappa3_hat_old = kappa3_hat;
    weighted_logllh_old = weighted_logllh;
    
        
    %% STEP 1: E-step
    for i = 1:k
        % test if there is a problem with the pdf calculation
        circular_pdf = circ_bvmpdf(data(:,1),data(:,2),...
            mu_hat(i), nu_hat(i), kappa1_hat(i),  kappa2_hat(i), kappa3_hat(i));
        if sum(isnan(circular_pdf)) || sum(isinf(circular_pdf))
            msg = ['Non feasible pdf for initial point ', ...
            ' parameters are: mu =  ',num2str(mu_hat(i)),', nu = ',...
            num2str(nu_hat(i)),', kappa1 = ',num2str(kappa1_hat(i)),....
            ', kappa2 = ',num2str(kappa2_hat(i)), ', kappa3 = ',num2str(kappa3_hat(i))];
            disp(msg);
            error('MATLAB:myCode:initPoint', msg);
        end
       
        posterior_probs(:,i) = prior_probs(i)*circ_bvmpdf(data(:,1),data(:,2),...
            mu_hat(i), nu_hat(i), kappa1_hat(i),  kappa2_hat(i), kappa3_hat(i));
    end
    if opts.noise
        posterior_probs(:,k+1) = prior_probs(k+1)*repmat(1/(2*pi)^2,numData,1);
    end
    % Normalize the posterior probs
    posterior_probs = posterior_probs./repmat(sum(posterior_probs,2),1,k + opts.noise);
    %[~,indx_membership] = max(posterior_probs,[],2);
    %% STEP 2: M-step
    for i = 1:k
        prior_probs(i) = mean(posterior_probs(:,i));
        % weighted sample by the posterior distribution
        [param_best, funcval_final, exitflag] = fmincon(@(x) ...
            - Loglikelihood(posterior_probs(:,i),data(:,1),data(:,2),x),...
            [mu_hat(i),nu_hat(i),kappa1_hat(i),kappa2_hat(i)],...%kappa3_hat(i)],...
            [],[],[],[],lb(i,:),ub(i,:),@confun,fmincon_opts);     
      
        % deterministic annealing
        kappa1_hat(i) = min(param_best(3),max_kappa) * (param_best(3) >= threshold_kappa) + ...
            min(max_kappa, max(param_best(3), kappa1_hat_old(i)*mult))*(param_best(3) < threshold_kappa);
        kappa2_hat(i) = min(param_best(4),max_kappa) * (param_best(4) >= threshold_kappa) + ...
            min(max_kappa, max(param_best(4), kappa2_hat_old(i)*mult))*(param_best(4) < threshold_kappa);
        mu_hat(i) = param_best(1);
        nu_hat(i) = param_best(2);
        %kappa1_hat(i) = param_best(3);
        %kappa2_hat(i) = param_best(4);
        kappa3_hat(i) = 0; %param_best(5);     
        weighted_logllh(i) = Loglikelihood(data(:,1),data(:,2),posterior_probs(:,i),...
            [mu_hat(i),nu_hat(i),kappa1_hat(i),kappa2_hat(i)]);
    end
    
    % rescale the uniform noise if it goes above 10%
    noise_threshold = 0.01;
    if opts.noise && sum(prior_probs(1:k)) < 1- noise_threshold
        prior_probs(1:k) = prior_probs(1:k)*(1-noise_threshold+noise_threshold*rand)/sum(prior_probs(1:k));
    end
    
    if opts.noise
        prior_probs(k+1) = 1 - sum(prior_probs(1:k));
    end
        
    %% Stopping criteria
        
    %% CHANGE THE LIKELIHOOD CONDITION
    %llh_change = norm(abs(weighted_logllh - weighted_logllh_old));
    llh_change = abs((sum(weighted_logllh) - sum(weighted_logllh_old))./sum(weighted_logllh_old));
    mu_change = norm(abs(mu_hat - mu_hat_old));
    nu_change = norm(abs(nu_hat - nu_hat_old));
    kappa1_change = norm(abs(kappa1_hat - kappa1_hat_old));
    kappa2_change = norm(abs(kappa2_hat - kappa2_hat_old));
    %kappa3_change = norm(abs(kappa3_hat - kappa3_hat_old));
    
    if llh_change  < opts.eps || (mu_change  < opts.eps && nu_change  < opts.eps ...
            && kappa1_change < opts.eps && kappa2_change  < opts.eps )%|| kappa3_change  < opts.eps
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

function [LLH] = Loglikelihood(phi, psi, pij, params)
    mu = params(1); nu = params(2); kappa1 = params(3); kappa2 = params(4); %kappa3 = params(5);
    % in case of independence kappa3 = 0        
    LLH = sum(pij) + sum(kappa1*cos(phi-mu)) + sum(kappa2*cos(psi-nu)) -...
      length(phi)*log(2*pi*besseli(0,kappa1)*2*pi*besseli(0,kappa2));     
end


% function [LLH] = Loglikelihood(phi, psi, pij, params)
%     mu = params(1); nu = params(2); kappa1 = params(3); kappa2 = params(4); kappa3 = params(5);
%     %H = kappa1*cos(phi-mu) + kappa2*cos(psi - nu) - kappa3*cos(phi-mu-psi+nu);
%     % in case of independence kappa3 = 0    
%     %fun = @(x, nu, kappa1, kappa2, kappa3) 2*pi*besseli(0,sqrt(kappa1.^2+kappa3.^2 ...
%     %-2*kappa1.*kappa3.*cos(x - nu))).*exp(kappa2.*cos(x-nu));
%     %Cc_inv = integral((@(x)fun(x, nu, kappa1, kappa2, kappa3)),0,2*pi);
%     %LLH = length(phi)*log(Cc_inv^-1) + sum(H) + sum(log(pij));
%     fun_kappa13 = @(ang) sqrt(kappa1.^2+kappa3.^2 -2*kappa1.*kappa3.*cos(ang - nu));
%     fun_kappa23 = @(ang) sqrt(kappa2.^2+kappa3.^2 -2*kappa2.*kappa3.*cos(ang - mu));
%     psi_nu = atan( - kappa3.*sin(psi - nu)./(kappa1 - kappa3.*cos(psi-nu)));
%     phi_mu = atan( - kappa3.*sin(phi - mu)./(kappa2 - kappa3.*cos(phi-mu)));
%     H = fun_kappa13(psi).*cos(phi - mu + psi_nu) + fun_kappa23(phi).*cos(psi - nu + phi_mu); 
%     if abs(kappa3) > 1e-5
%         LLH = sum(pij) + sum(H) - length(phi)*(log(2^2*pi^2)) - ...
%         sum(log(besseli(0,fun_kappa13(phi)))) - sum(log(besseli(0,fun_kappa23(phi)))); 
%     else
%         LLH = sum(pij) + sum(kappa1*cos(phi-mu)) + sum(kappa2*cos(psi-nu)) -...
%             length(phi)*(2^2*pi^2*besseli(0,kappa1)*besseli(0,kappa2));
%     end 
% end

function [c, ceq] = confun(params)
kappa1 = params(3); kappa2 = params(4); %kappa3 = params(5);
c = [];%skappa3 - kappa1*kappa2/(kappa1+kappa2);
ceq = [];
end
