%% function [W] = buildW_pmi(f_maps,rf,p,opts,samples)
% builds an affinity matrix W for image im based on PMI
% 
% INPUTS
%  f_maps   - NxMxF array of F feature maps for an NxM image
%  rf       - the learned random forest for approximating PMI (unused if ~opts.approximate_PMI)
%  p        - P(A,B) (unused if opts.approximate_PMI)
%  mixture_params - P(A,B) from mixture model (unused if
%  opts.approximate_PMI)
%  opts     - parameter settings (see setEnvironment)
%  samples   - either the number of samples from the full affinity matrix to
%               compute, or the indices of the full affinity matrix to compute, or empty,
%               in which case the full affinity matrix is computed
%
% OUTPUTS
%  W - affinity matrix
%
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2014 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------
% Luong Nguyen 10/06/14: change line 32, add rad,rad_inner
function [W] = buildW_pmi(f_maps,rf,p,mixture_params, which_feature, opts,samples)
    
    if (~exist('samples','var'))
        samples = [];
    end
    
    im_size = size(f_maps(:,:,1));
    
    %% get local pixel pairs
    if (isempty(samples) || size(samples,2)==1)
        [ii,jj] = getLocalPairs(im_size,opts.localPairs.rad,opts.localPairs.rad_inner,samples);
    else
        ii = samples(:,1);
        jj = samples(:,2);
    end
    
    %% initialize affinity matrix
    Npixels = prod(im_size);
    W = sparse(double(ii),double(jj),0,Npixels,Npixels);
    
    %% extract features F
    [F,F_unary] = extractF(f_maps,ii,jj,opts);
    
    %% evaluate affinities
    if (opts.approximate_PMI)
        w = exp(fastRFreg_predict(F,rf));
        %w = fastRFreg_predict(F,rf);
    else
        if strcmp(which_feature,'hue opp')
            pmi = evalPMI_theta(F,mixture_params,opts);
        else
            pmi = evalPMI(p,F,F_unary,ii,jj,opts);
        end
        w = exp(pmi);
        %w = pmi;
    end
    
    %%
    W2 = sparse(double(ii),double(jj),w,Npixels,Npixels);
    W = W+W2;
    W = (W+W'); % we only computed one half of the affinities, now assume they are symmetric
end