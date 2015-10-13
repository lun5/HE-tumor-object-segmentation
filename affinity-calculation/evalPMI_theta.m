% calculating the PMI for H&E hue
% Luong Nguyen 5/15/2015
% INPUT: h&e hue values stored in F(:,1) and F(:,2)
%        mixture_params: paramesters of mixture of bivariate VM 
%        containing: mu,nu, kappa1, kappa2, kappa3, prior_probs
%       opts: from setEnvironment_affinity
function [pmi,pJoint,pProd] = evalPMI_theta(F,mixture_params,opts)
   %% evaluate p(A,B)
    % parameters of the mixture model
    params = mixture_params.params;
    prior_probs = mixture_params.prior_probs;
    init_params = mixture_params.init_params;
    mu = params.mu; nu = params.nu; 
    %kappa1 = params.kappa1; kappa2 = params.kappa2; kappa3 = params.kappa3;
 
%% evaluate these joint distribution at the sampled points
    prc = 5;
    % cap the joint distribution
%     mult = 1;%1.5;
%     pJoint_max = mult.*max(jointDist(mu(1),nu(1),params,prior_probs),...
%         jointDist(mu(2),nu(2),params,prior_probs));
%     pJoint = min(pJoint,pJoint_max);
    ratio_white = min(jointDist(mu(1),nu(1),params,prior_probs),...
         jointDist(mu(2),nu(2),params,prior_probs))./(jointDist(mu(3),nu(3),params,prior_probs));
    ratio_white = min(ratio_white,1);
    prior_probs(3) = prior_probs(3)*ratio_white;    
    pJoint = jointDist(F(:,1), F(:,2), params, prior_probs);
    if (opts.model_half_space_only)
        pJoint = pJoint./2; % divided by 2 since we only modeled half the space
    end
    %% evaluate p(A)p(B)
    %pMarg_phi = marginalDist(F(:,1), params, prior_probs, 1);
    % cap the marginal distribution
%     mult = 1;
%     pMarg_max = mult.*max(marginalDist(mu(1),init_params),marginalDist(mu(2),init_params));
%     pMarg_phi = min(pMarg_phi,pMarg_max);pMarg_psi = min(pMarg_psi,pMarg_max);
    ratio_white = min(marginalDist(init_params.theta_hat(1),init_params),...
        marginalDist(init_params.theta_hat(2),init_params))./...
        (1.5*marginalDist(init_params.theta_hat(3),init_params));
    ratio_white = min(ratio_white,1);
    init_params.prior_probs(3) = ratio_white * init_params.prior_probs(3);
    pMarg_phi = marginalDist(F(:,1), init_params);
    pMarg_psi = marginalDist(F(:,2), init_params);

    pProd = pMarg_phi.*pMarg_psi;
    reg = prctile(nonzeros(pProd),prc);
    pProd = pProd+reg; pJoint = pJoint + reg;

    %% calculate pmi
    %pmi = ((pJoint+reg).^(opts.joint_exponent))./pProd;
    pmi = log(pJoint)*opts.joint_exponent - log(pProd);
end

% a function for the joint distribution
function pJoint = jointDist(phi, psi, params, prior_probs)
    mu = params.mu; nu = params.nu; kappa1 = params.kappa1; 
    kappa2 = params.kappa2; kappa3 = params.kappa3;
    numClusters = length(mu);
    numPoints = length(phi);
    pJoint = zeros(numPoints,1);
    for i = 1:numClusters       
        pJoint = pJoint + prior_probs(i)*circ_bvmpdf(phi,psi,mu(i),nu(i),kappa1(i),kappa2(i),kappa3(i));
        %pJoint = Cc_inv(i).^-1 * exp(kappa1(i)*cos(phi-mu(i)) +...
        %    kappa2(i)*cos(psi-nu(i)) - kappa3(i)*cos(phi-mu(i) -psi+nu(i)));
    end
end

function pMarginal = marginalDist(ang,init_params)
    mu_hat = init_params.theta_hat;
    kappa_hat = init_params.kappa_hat;
    prior_probs = init_params.prior_probs;
    numCluster = length(mu_hat);
    pMarginal = zeros(size(ang));
    for i = 1:numCluster
        pMarginal = pMarginal + prior_probs(i)*circ_vmpdf(ang,mu_hat(i),kappa_hat(i));
    end
end

% function pMarginal = marginalDist(ang, params, prior_probs,dim)
%     mu = params.mu; nu = params.nu; kappa1 = params.kappa1; 
%     kappa2 = params.kappa2; kappa3 = params.kappa3;
%     numPoints = length(ang); numClusters = length(mu);
%     pMarginal = zeros(size(ang));
%     %if size(ang,1) < size(ang,2); ang = ang'; end
%     fun_kappa13 = @(x, kappa1, kappa2, kappa3, nu) ...
%         sqrt(kappa1.^2+kappa3.^2 -2*kappa1.*kappa3.*cos(x - nu));
%     fun_kappa23 = @(x, kappa1, kappa2, kappa3, mu) ...
%         sqrt(kappa2.^2+kappa3.^2 -2*kappa2.*kappa3.*cos(x - mu));
%     
%      % normalizing factors
%     fun_Cc_inv_psi = @(x, nu, kappa1, kappa2, kappa3) 2*pi*besseli(0,...
%         sqrt(fun_kappa13(x,kappa1,kappa2,kappa3,nu))).*exp(kappa2.*cos(x-nu));
%     fun_Cc_inv_phi = @(x, mu, kappa1, kappa2, kappa3) 2*pi*besseli(0,...
%         sqrt(fun_kappa23(x,kappa1,kappa2,kappa3,mu))).*exp(kappa1.*cos(x-mu));
% 
%     if dim == 1
%         for i = 1:numClusters
%             Cc_inv = integral((@(x)fun_Cc_inv_phi(x, mu(i), kappa1(i), kappa2(i), kappa3(i))),0,2*pi);
%             pMarginal = pMarginal + prior_probs(i).*Cc_inv^(-1)*2*pi.*besseli(0,...
%             sqrt(fun_kappa23(ang,kappa1(i),kappa2(i), kappa3(i), mu(i))))...
%             .*exp(kappa1(i).*cos(ang - repmat(mu(i),size(ang))));
%         end
%     elseif dim == 2
%         for i = 1:numClusters
%             Cc_inv = integral((@(x)fun_Cc_inv_psi(x, nu(i), kappa1(i), kappa2(i), kappa3(i))),0,2*pi);
%             pMarginal = pMarginal + prior_probs(i)*Cc_inv^(-1)*2*pi.*besseli(0,...
%             sqrt(fun_kappa13(ang,kappa1(i),kappa1(i), kappa3(i), nu(i))))...
%             .*exp(kappa2(i)*cos(ang - repmat(nu(i),size(ang))));
%         end
%     else
%         
%     end
% end

