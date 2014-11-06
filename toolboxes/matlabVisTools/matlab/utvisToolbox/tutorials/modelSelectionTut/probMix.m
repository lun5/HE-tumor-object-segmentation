function [pSum, p] = probMix(model, xRange)
%  [pSum, p] = probMixComp(model, xRange)
%
% Compute the probability of each component of a given mixture model
%  to give a response within a particular range of the independent
%  variable x in xRange = [xMin, xMax].
%  The probability that mixture component k lands in this interval is returned
%  in p(k).  pSum = sum(p);

  h = xRange(2) - xRange(1);
  pSum = 0.0;
  nComp = model.nInliers+1;
  p = zeros(nComp,1);
  sqrtHalf = sqrt(0.5);

  %%% Outlier component 
  p(nComp) = h * model.outLike * model.mix(nComp);
  pSum = p(nComp);

  %%% Gaussian model components 
  for  m=1:model.nInliers
    if (model.sig(m) > 2.0 * h)
      %% Use a trapezoidal approx to the integral 
  p(m) = (h/2.0) * model.mix(m) * ...
	( gauss1D(xRange(1), model.mean(m), model.sig(m)) ...
	  + gauss1D(xRange(2), model.mean(m), model.sig(m)) );
	 
    else
      % When sig < 2h the trapezoidal approx will be
      % poor.  Use a refined estimate of the integral provided
      % by the built-in error function.
      dev0 = (xRange(1) - model.mean(m))/(model.sig(m));
      dev1 = (xRange(2) - model.mean(m))/(model.sig(m));
      if (abs(dev0) < 6.0 | abs(dev1) < 6.0 | dev1 * dev0 < 0) 
	p(m) = model.mix(m) * ...
	  (erf(dev1 * sqrtHalf) - erf(dev0 * sqrtHalf))/2.0;
      else
	p(m) = 0.0;
      end
    end
    pSum = pSum + p(m);
  end %% loop over Gaussian mixture components 
  return;

