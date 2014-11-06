function  p = priorMixModel(model) 
% p = priorMixModel(model)
%
% Use a uniform prior over mixing coeffs, mean, and log(sigma).
% Return the probability, according to this prior, of selecting
% parameters from a particular range of values.  The range
% is selected to be roughly the resolution of the fitted models.
% That is, the region over which the total log likelihood is
% roughly constant.
%
% Input:
%  model - mix model struct
% Output:
%  p - prior probability

  p = 1.0;
  for m=1:model.nInliers
    %%% resolution of mixing component is +-0.10 / 1.0
    p = p * 2*0.1;  
    %%% sigma = e^s in [1 64],
    %   with s uniform in [log(1.0), log(64.0)]. 
    %   Also sigma is resolved to a factor of 2. 
    %   Therfore, delta_s is resolved to +-log(2.0) from a uniform
    %   distribution  over [0, log(64)]
    p = p * 2.0 * log(2.0)/log(64.0); 
    % mean uniform in [0 256], mean resolved to +- sigma */
    p = p * 2.0 * ( model.sig(m)) / 256.0;
  end 
  return;

