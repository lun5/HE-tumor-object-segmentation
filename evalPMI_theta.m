% calculating the PMI
% INPUT: joint distribution parameters
%       save_struct: a structure with the following fields
%           params.mu,nu, kappa1, kappa2, kappa3
%           prior_probs
%       opts.rho: power of joint distribution
% F_unary,A_idx,B_idx,
function [pmi,pJoint,pProd] = evalPMI_theta(F,mixture_params,opts)
   %% evaluate p(A,B)
    reg = opts.p_reg;
    %tol = opts.kde.kdtree_tol;
    
    % parameters of the mixture model
    params = mixture_params.params;
    prior_probs = mixture_params.prior_probs;
    mu = params.mu; nu = params.nu; kappa1 = params.kappa1; 
    kappa2 = params.kappa2; kappa3 = params.kappa3;
    % normalizing factors
    fun_Cc_inv = @(x, nu, kappa1, kappa2, kappa3) 2*pi*besseli(0,sqrt(kappa1.^2+kappa3.^2 ...
    -2*kappa1.*kappa3.*cos(x - nu))).*exp(kappa2.*cos(x-nu));
    numClusters = length(mu); Cc_inv = zeros(size(mu));
    for i = 1:numClusters
        Cc_inv(i) = integral((@(x)fun_Cc_inv(x, nu(i), kappa1(i), kappa2(i), kappa3(i))),0,2*pi);
    end
        
    % evaluate these joint distribution at the sampled points
    pJoint = reg + jointDist(F(:,1), F(:,2), params, prior_probs)/2; % divided by 2 since we only modeled half the space

    %% evaluate p(A)p(B)
    pMarg_phi = marginalDist(F(:,1), params, prior_probs,Cc_inv, 1);
    pMarg_psi = marginalDist(F(:,2), params, prior_probs,CC_inv, 2);
    pProd = pMarg_phi.*pMarg_psi+reg;

    %% calculate pmi
    pmi = log((pJoint.^(opts.joint_exponent))./pProd);
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
    end
end

% a function for the marginal distribution
% dim: dimension along which we marginalize
% dim = 1: pMarginal gives mariginal distribution of phi
% dim = 2: pMarginal gives mariginal distribution of psi
function pMarginal = marginalDist(ang, params, prior_probs, Cc_inv dim)
    mu = params.mu; nu = params.nu; kappa1 = params.kappa1; 
    kappa2 = params.kappa2; kappa3 = params.kappa3;
    numPoints = length(ang); numClusters = length(mu);
    pMarginal = zeros(numPoints,1);

    fun_kappa13 = @(x, kappa1, kappa2, kappa3, nu) ...
        sqrt(kappa1.^2+kappa3.^2 -2*kappa1.*kappa3.*cos(x - nu));
    fun_kappa23 = @(x, kappa1, kappa2, kappa3, mu) ...
        sqrt(kappa2.^2+kappa3.^2 -2*kappa2.*kappa3.*cos(x - mu));
    if dim == 1
        for i = 1:numClusters
            pMarginal = pMarginal + prior_probs(i)*Cc_inv(i)^(-1)*2*pi*besseli(0,...
            sqrt(fun_kappa23(ang,kappa1(i),kappa2(i), kappa3(i), mu(i))))...
            .*exp(kappa1(i).*cos(ang - repmat(mu(i),size(ang))));
        end
    elseif dim == 2
        for i = 1:numClusters
            pMarginal = pMarginal + prior_probs(i)*Cc_inv(i)^(-1)*2*pi*besseli(0,...
            sqrt(fun_kappa13(ang,kappa1(i),kappa1(i), kappa3(i), nu(i))))...
            .*exp(kappa2(i)*cos(ang - repmat(nu(i),size(ang))));
        end
    else
        
    end
end

