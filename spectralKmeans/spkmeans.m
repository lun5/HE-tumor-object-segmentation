%% [idxbest, mu_hat_cart] = spkmeans(X_cart,k,opts)% 
% INPUTS ONLY IN 2D!!!!!
%  X             - data in cartesian space of size mxd. d: dimension, m:
%                - number of data points. Will convert to polar
%                    
%  k             - number of clusters to find.
%  opts          - parameter settings (maxiter, stopping criteria)
%
% OUTPUTS
%  mu_hat_cart - estimated centroid cluster in Cartesian space
%  kappa_hat     - estimated concentration parameters
%  mu_hat_cart   - estimated means in cartesian space
%  posterior_probs - component membership P(in kth component| data i)
%
% -------------------------------------------------------------------------
% HE segmentation toolbox
% Luong Nguyen, 2014 [lun5@pitt.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

function [idxbest, centroids_cart] = spkmeans(X_cart,k,opts)
    
    opts_default.maxiter = 1000;
    opts_default.eps = 1e-3; % threshold for likelihood convergence
    %opts_default.minNumElt = size(X_cart,2)/(k+100);
    if nargin <3
        opts = opts_default;
    elseif nargin <2
        error('Function needs at least 2 inputs: data, number of components');
    end
    
    if ~exist('opts','var')
        opts = opts_default;
    end
    
    if ~exist('opts.maxiter','var')
        opts.maxiter = opts_default.maxiter;
    end
    
    if ~exist('opts.eps','var');
        opts.eps = opts_default.eps;
    end
% 
%     if ~exist('opts.minNumElt','var');
%         opts.minNumElt = opts_default.minNumElt;
%     end

    % Set 'm' to the number of data points.
    numData = size(X_cart, 1);
    % set 'd' to the dimension
    d = size(X_cart,2); % for now it's just 2D
    %X_polar = atan2(X_cart(:,2),X_cart(:,1));

    %% STEP 1: Initialization
    centroids_cart = zeros(d,k); % centroids in cartesian
    mu_hat_polar = zeros(1,k); % centroids in polar
    
    %% randomly assign mu
    for i = 1:k
        mu_hat_polar(i) = (i-1)*pi/k;
        centroids_cart(:,i) = [cos(mu_hat_polar(i)); sin(mu_hat_polar(i))];
    end

    % indx for the clusters
    idxbest = zeros(numData,1);
    
for iter = 1: opts.maxiter
    
    mu_hat_cart_old = centroids_cart; 
    
    %% STEP 1: E-step
    distance_to_centroids = X_cart*centroids_cart;
    [~,idxbest] =  max(distance_to_centroids,[],2);
    %% STEP 2: M-step
    for j = 1:k
        indx_cluster = idxbest == j;
        centroids_cart(:,j) = sum(X_cart(indx_cluster,:),1)'/norm(sum(X_cart(indx_cluster,:),1));
    end
        
    %% Stopping criteria
    centroid_change = centroids_cart - mu_hat_cart_old;
    if norm(centroid_change(:)) < opts.eps
        break;
    end
    % stop if the number of a 
    
    
end

if iter == opts.maxiter
    sprintf('The algorithm does not converge after maxiter = %d',opts.maxiter);
end
  
end