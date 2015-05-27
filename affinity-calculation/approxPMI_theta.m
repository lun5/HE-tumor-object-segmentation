% calculating the PMI
% INPUT: joint distribution parameters
%       save_struct: a structure with the following fields
%           params.mu,nu, kappa1, kappa2, kappa3
%           prior_probs
%       opts.rho: power of joint distribution
% F_unary,A_idx,B_idx,
function [pmi,pJoint,pProd] = approxPMI_theta(F,mixture_params,opts)
   %% evaluate p(A,B)
    %reg = opts.p_reg;
    %tol = opts.kde.kdtree_tol;
    
    % parameters of the mixture model
    params = mixture_params.params;
    prior_probs = mixture_params.prior_probs;
    mu = params.mu; nu = params.nu; kappa1 = params.kappa1; 
    kappa2 = params.kappa2; kappa3 = params.kappa3;
    %% normalizing factors
    fun_Cc_inv = @(x, nu, kappa1, kappa2, kappa3) 2*pi*besseli(0,sqrt(kappa1.^2+kappa3.^2 ...
    -2*kappa1.*kappa3.*cos(x - nu))).*exp(kappa2.*cos(x-nu));
    numClusters = length(mu); numPoints = size(F,1);
    Cc_inv = zeros(size(mu));
    for i = 1:numClusters
       if i > 6
           Cc_inv(i) = Cc_inv(i-3); continue
       end
       Cc_inv(i) = integral((@(x)fun_Cc_inv(x, nu(i), kappa1(i), kappa2(i), kappa3(i))),0,2*pi);       
    end
    % evaluate these joint distribution at the sampled points
    %prc = 5;
    logpJoint = zeros(numPoints,1);
    
    for i = 1:numClusters
        logpJoint = logpJoint + prior_probs(i).*(-log(Cc_inv(i)) + ...
            kappa1(i).*cos(F(:,1) - mu(i)) + kappa2(i).*cos(F(:,2)- nu(i)) - ...
            kappa3(i).*cos(F(:,1) - mu(i) - F(:,2) + nu(i))); 
    end
    pJoint = exp(logpJoint);
    if (opts.model_half_space_only)
        pJoint = pJoint./2; % divided by 2 since we only modeled half the space
        logpJoint = log(pJoint+eps);
    end

    %% evaluate p(A)p(B)
    logPMarg_phi = zeros(numPoints,1); logPMarg_psi = zeros(numPoints,1); 
    fun_kappa13 = @(x, kappa1, kappa2, kappa3, nu) ...
        sqrt(kappa1.^2+kappa3.^2 -2*kappa1.*kappa3.*cos(x - nu));
    fun_kappa23 = @(x, kappa1, kappa2, kappa3, mu) ...
        sqrt(kappa2.^2+kappa3.^2 -2*kappa2.*kappa3.*cos(x - mu));

    for i = 1:numClusters
        logPMarg_phi = logPMarg_phi + prior_probs(i)*(-log(Cc_inv(i)) + log(2*pi) + ...
            + log(besseli(0,sqrt(fun_kappa23(F(:,1),kappa1(i),kappa2(i), kappa3(i), mu(i))))) + ...
            + kappa1(i)*cos(F(:,1) - mu(i)));
        logPMarg_psi = logPMarg_psi + prior_probs(i)*(-log(Cc_inv(i)) + log(2*pi) + ...
            + log(besseli(0,sqrt(fun_kappa13(F(:,2),kappa1(i),kappa2(i), kappa3(i), nu(i))))) + ...
            + kappa2(i)*cos(F(:,2) - nu(i)));
    end
    
    pProd = exp(logPMarg_phi + logPMarg_psi);
    %% calculate pmi
    pmi = logpJoint*opts.joint_exponent - logPMarg_phi - logPMarg_psi;
end
