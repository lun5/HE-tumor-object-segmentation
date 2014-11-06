function [] = plotLinBCC(C, vm, nPlot)
%  plotLinBCC(C, vm, nPlot): Plot brightness constancy constraints
%
%  Input:
%     C(1:3, k) contains motion constraint data (gradI , deltaI)^T
%     vm  plot velocity within an axis of range +-vm in both vx and vy
%     nPlot randomly sample at most nPlot of the constraints.
  
  if nargin < 2
    vm = 2.0;
  end
  if nargin < 3
    nPlot = 200;
  end
  h = gcf;
  axis([-vm vm -vm vm]); axis 'square'; 
  set(gca,'Ydir','reverse'); 
  cropBox = [-vm -vm vm vm];
  n = size(C,2);
  idx = randperm(n);
  idx = idx(1:min(length(idx), nPlot));
  
  for k = 1:length(idx)
    endPnts = cropLineInBox(C(1:2,idx(k)), C(3,idx(k)), cropBox);
    if ~isnan(endPnts(1))
      set(line(endPnts(:,1), endPnts(:,2)), 'Color', ...
                        [0 0 0]);
    end
  end

  return;