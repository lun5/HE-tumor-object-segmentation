function [totalLogLike, own, dataLike] = dataLogLikeMixModel(data, model)
% [totalLogLike, own, dataLike] = dataLogLikeMixModel(data, model)
%
% Return the log likelihood of a 1D data set given a mixture
% model.  Also, own is a (length(data), model.nInliers+1) matrix
% containing the ownership probabilities for each mixture model
% component and each data item.  dataLike is the likelihood
% of each data item separately (totalLogLike=sum(log(dataLike))

% constants
sqrt2pi = sqrt(2*pi);
N = length(data(:));

if (size(model.mix,2) > 1)
  model.mix = model.mix';
end
layers = model.nInliers+1;


totalLogLike = 0.0;
own = zeros(N,layers);

  % E- step
  % Compute component distribution likelihoods
  compLike = zeros(length(data), layers);
  for k = 1:model.nInliers
    tmp = (data - model.mean(k))/model.sig(k);
    tmp = -0.5 * (tmp .^ 2);
    compLike(:,k) = exp(tmp)/(sqrt2pi*model.sig(k));
  end
  compLike(:, layers) = model.outLike*ones(length(data),1);

  % Data likelihood is the sum of the mixing probs times the
  % component likelihoods.
  dataLike = compLike * model.mix;
  
  % Total (log) likelihood of data set
  totalLogLike = sum(log(dataLike));
  
  % ownership probabilities
  for k = 1:layers
    own(:,k) = (model.mix(k)*compLike(:,k))./dataLike;
  end

return












