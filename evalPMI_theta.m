% calculating the PMI
% INPUT: joint distribution parameters
%       save_struct: a structure with the following fields
%           params.mu,nu, kappa1, kappa2, kappa3
%           prior_probs
%       opts.rho: power of joint distribution

function [pmi,pJoint,pProd] = evalPMI_theta(mixture_params,F,F_unary,A_idx,B_idx,opts)
   %% evaluate p(A,B)
    reg = opts.p_reg;
    %tol = opts.kde.kdtree_tol;
    
    % parameters of the mixture model
    params = mixture_params.params;
    prior_probs = mixture_params.prior_probs;
    
    % reflect all of these to have the full space of the joint distribution
    
    % a function for the joint distribution 
    
    % a function for the marginal distribution
    
    % evaluate these functions at the samples
    pJoint = reg + evaluate_batches(p,F',tol)/2; % divided by 2 since we only modeled half the space

    %% evaluate p(A)p(B)
    N = floor(size(F,2)/2); assert((round(N)-N)==0);
    p2_1 = marginal(p,1:N);
    p2_2 = marginal(p,N+1:(2*N));
    p2 = joinTrees(p2_1,p2_2,0.5);
    pMarg = zeros(size(F_unary,1),1);
    ii = find(~isnan(F_unary(:,1))); % only evaluate where not nan (A_idx and B_idx will only refer to non-nan entries)
    pMarg(ii) = evaluate_batches(p2,F_unary(ii,:)',tol);
    pProd = pMarg(A_idx).*pMarg(B_idx)+reg;

    %% calculate pmi
    pmi = log((pJoint.^(opts.joint_exponent))./pProd);

end