function [fit, own, totalLogLike, its] = fitMixModel(data, guess,...
                    emOpt, minSigma)
%[fit, own, totalLogLike, its] = fitMixModel(data, guess, options, minSigma)
% 
% Fit Gaussian + outlier mixture model to 1D (gray level) data.
% Input:
%   data  - N x 1 array of data 
%   guess - mixture model struct for initial guess, see mixScript.m
%   emOpt:  
%     emOpt.tolLike  the step tolerance in the overall data likelihood
%                    for the stopping criterion.
%     emOpt.maxIts   the maximum number of iterations.
%   minSigma - minimum value of component distribution standard deviations.
% Output
%   fit - the fitted mixture model struct
%   own - N x (guess.nInliers + 1) data ownership matrix
%   totalLogLike - log likelihood of data according to fit model (a scalar)
%   its - number of EM iterations used to obtain the fit model
%--------------------------------------------------------------
% CC, ADJ - Nov, 01

% Optional input default values
if (nargin<4)
  minSigma = 1;
end

test = 0;
if emOpt.tolLike > 0.0
  test = 1;	% Test log likelihood for termination
end

iter = emOpt.maxIts;

% constants
sqrt2pi = sqrt(2*pi);
N = length(data(:));

% Copy initial guess to fit model.
fit = guess;

if length(data) == 0
  return
end

% Do EM iterations.
layers = fit.nInliers+1;
totalLogLike = zeros(iter, 1);
for its = 1:iter

  % E- step
  % Compute component distribution likelihoods
  compLike = zeros(length(data), layers);
  for k = 1:fit.nInliers
    tmp = (data - fit.mean(k))/fit.sig(k);
    tmp = -0.5 * (tmp .^ 2);
    compLike(:,k) = exp(tmp)/(sqrt2pi*fit.sig(k));
  end
  compLike(:, layers) = fit.outLike*ones(length(data),1);

  % Data likelihood is the sum of the mixing probs times the
  % component likelihoods.

  dataLike = compLike * ((fit.mix)');
  
  % Total (log) likelihood of data set
  totalLogLike(its) = sum(log(dataLike));
  
  % ownership probabilities
  for k = 1:layers
    own(:,k) = (fit.mix(k)*compLike(:,k))./dataLike;
  end

  % M-step

  % Update mixing coeffs m
  dataOwn = sum(own,1)';
  fit.mix = dataOwn'/sum(dataOwn);
  
  % Update means
  for k = 1:layers-1,
    if dataOwn(k) > 0.0
      fit.mean(k) = sum(own(:,k).*data)/dataOwn(k);
    end
  end
  
  % update sigma
  for k = 1:layers-1,
    if dataOwn(k) > 0.0
      p = fit.mean(k);
      fit.sig(k) = sqrt( sum(own(:,k).*((data-p).^2))/dataOwn(k) );
    end
    fit.sig(fit.sig<minSigma) = minSigma;
  end

  % check if log-likelihood increased by significant amount
  if test
    if (its > 1 && abs(totalLogLike(its) - totalLogLike(its-1)) < emOpt.tolLike)
      break;
    end
  end
  
end
return












