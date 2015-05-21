function plotData(binHist, dataHist, yMax)
% plotData(binHist, dataHist, yMax)
%
% Useage:
%   See mixScript.m
% Plot histograms for samples from a 1D gray-level mix model.
% Optional:
%   yMax: Current maximum for y scaling

  if (nargin < 2) 
    disp('plotData requires >=2 argument(s)');
    return;
  end;
  if (nargin < 3) 
     yMax = 0.0; 
  end;

  % Choose upper bound for y-axis in the plots
  yMax = max([max(dataHist), yMax]);

  if yMax > 0.1
    yMax = ceil(yMax/0.05) * 0.05;
  else
    yMax = ceil(yMax/0.02) * 0.02;
  end

  bar(binHist,dataHist(:,1),1,'g'); axis([0 255 0 yMax]), axis square;
  xlabel('Grey level'); ylabel('Frequency');
  title('Data Samples');

