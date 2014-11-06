function x = plotMix(modelHist, yMax)
% PLOTMIX: Plot histograms for a 1D gray-level mix model
% 
% Useage:
%   See mixScript.m
% Input:
%   modelHist: model histogram as defined in mixScript.m
%   yMax (optional 0.0) minimum y-axis height to use.
% Returns x positions plotted

  if nargin < 2
    yMax = 0.0;
  end

  nHist = size(modelHist, 1);
  dx=255/(nHist-1);
  x = 0:dx:255;
  % Choose upper bound for y-axis in the plots
  yMax = max([max(modelHist) yMax]);
  if yMax > 0.1
    yMax = ceil(yMax/0.05) * 0.05;
  else
    yMax = ceil(yMax/0.02) * 0.02;
  end

  plot(x, modelHist, 'LineWidth', 2.0), axis([0 255 0 yMax]), axis square;
  xlabel('Grey level'); ylabel('Probability');
  title('Mixture Model');
  
