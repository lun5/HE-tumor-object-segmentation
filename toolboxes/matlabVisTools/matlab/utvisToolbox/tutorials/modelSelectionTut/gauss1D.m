function p = gauss1D(x, mn, sig)
%  p = gauss1D(x, mn, sig)
%
% Evaluate a 1D Normal density function at a vector
% of positions x.  The mean and variance of the Normal
% are given by mn, and sig, respectively.

  scl = sqrt(2.0 * pi)*sig;
  dev = (x - mn)/sig;
  dev = -0.5 * (dev .* dev);
  p = exp(dev)/scl; 
