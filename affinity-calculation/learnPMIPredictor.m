%% function [rf] = learnPMIPredictor(f_maps,p,opts)
% learns a random forest rf that approximates PMI_{\rho}(A,B)
% 
% INPUTS
%  f_maps   - NxMxF array of F feature maps for an NxM image
%  p        - model for P(A,B) (Eqn. 1 in paper)
%  opts     - parameter settings (see setEnvironment)
%
% OUTPUTS
%  rf       - random forest that approximates PMI_{\rho}(A,B)
%
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2014 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------
% Luong Nguyen 10/06/14: change line 32, add rad,rad_inner
% Luong Nguyen 5/7/15: build prediction model for hue opponent color
function [rf] = learnPMIPredictor(f_maps,p,mixture_params, which_feature, opts)
    %%
    Nsamples = opts.PMI_predictor.Nsamples_learning_PMI_predictor;
    im_size = size(f_maps(:,:,1));
    
    %% get all local pairs
    [ii,jj] = getLocalPairs(im_size,opts.localPairs.rad,opts.localPairs.rad_inner,Nsamples);
    
    %% extract features
    [F,F_unary] = extractF(f_maps,ii,jj,opts);
    
    %% check input
    if isempty(mixture_params) && strcmp(which_feature,'hue opp')
        error('Need mixture model parameters for hue opp');
    elseif isempty(p) && ~strcmp(which_feature,'hue opp')
        error('Need kde tree for non hue opponent color feature');     
    end
    
    %% evaluate affinities based on PMI
    if strcmp(which_feature,'hue opp')
        [pmi,pJoint,~] = evalPMI_theta(F,mixture_params,opts);
    else
        [pmi,pJoint,~] = evalPMI(p,F,F_unary,ii,jj,opts);
    end    
    %% learn PMI predictor: g(F) --> PMI
    Ntrees = opts.PMI_predictor.Ntrees;
    if strcmp(opts.affinityFunction,'PMI')
        rf = fastRFreg_train(F,pmi,Ntrees);
    elseif strcmp(opts.affinityFunction,'PJoint')
        rf = fastRFreg_train(F,pJoint,Ntrees);
    end
end