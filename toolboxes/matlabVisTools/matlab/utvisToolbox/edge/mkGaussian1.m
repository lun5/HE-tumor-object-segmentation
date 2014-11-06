% IM = mkGaussian1(SIZE, VARIANCE, MEAN, AMPLITUDE)
% 
% Compute a matrix with dimension SIZE (a scalar)
% containing a Gaussian function, centered at pixel position
% specified by MEAN (default = (size+1)/2), with given VARIANCE (a scalar).
% Default = (min(size)/6)^2)
% AMPLITUDE='norm' (default) ensures the sum of the weights is 1.

% modified from mkGaussian: Eero Simoncelli, 6/96.
% by ADJ, 10/01.


function [res] = mkGaussian(sz, var, mn, ampl)

sz = sz(:);
if (size(sz,1) ~= 1)
  sz = sz(1);
end

%------------------------------------------------------------
%% OPTIONAL ARGS:

if (exist('var') ~= 1)
  var = (sz(1)/6)^2;
end

if (exist('mn') ~= 1)
  mn = (sz+1)/2;
end

if (exist('ampl') ~= 1)
  ampl = 'norm';
end

%------------------------------------------------------------

xramp = [1:sz(1)] - mn;
res = exp(xramp.^2/(-2 * var(1)));

if (strcmp(ampl, 'norm'))
  res = res/sum(res);
end
