function data = sampleMix(model, nSamp)
%  data = sampleMix(model, nSamp)
%
% Compute nSamp independent samples from the 1D Gaussian mixture
% model specified by the struct model.  Assume that the outliers range
% over 0 to 1/outLike.

  cummMix = zeros(size(model.mix));
  for k = 1:length(cummMix)
    if k == 1
      cummMix(1) = model.mix(1);
    else
      cummMix(k) = model.mix(k) + cummMix(k-1);
    end
  end

  % Choose the component distribution
  u = rand(nSamp,1);
  comp = ones(nSamp, 1);
  for k = 1:model.nInliers
    indices = u>cummMix(k);
    comp(indices) = comp(indices) + 1;
  end
  
  data = zeros(nSamp,1);
  for c = 1:model.nInliers
    indices = (comp == c);
    data(indices) = model.sig(c) * randn(sum(indices),1) + model.mean(c);
  end
  
  c = model.nInliers + 1;
  indices = (comp == c);
  data(indices) = rand(sum(indices),1) / model.outLike;
  return;
  
