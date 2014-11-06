function [x, k] = prctile(d, pct)
  %% return point x and index k close to top pct percentile of 1D data d
  q = sort(d);
  k = floor(length(d) * (pct/100));
  k = max(k,1);
  x = q(k);