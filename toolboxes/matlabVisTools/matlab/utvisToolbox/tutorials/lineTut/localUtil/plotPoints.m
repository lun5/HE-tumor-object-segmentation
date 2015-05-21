function [plotBox] = plotPoints(p, h, lbl, targetLbl, plotBox)
% PLOTPOINTS Plot labelled image points, and color them
% according according to their labels.
% The points are plotted with x to the right, y down.
% (This is slow, don't use for too many thousands of points.)
% Arguments: p(:, 2) = (x1 y1; x2 y2; ...) point locations
% Optional Arguments:
%   h figure handle, a new figure is created if absent.
%   lbl(:)  If lbl is a vector then it provides the 
%           labels (integers>=1) for each of the points, i.e. 
%           lbl(i) is the label for point p(i,:)
%           The point is plotted in the color determined
%           by the lbl(i) entry of the current colormap,
%           with large values of lbl(i) wrapped into the range
%           of the current colormap.
%           If lbl is a scalar, then all points are plotted
%           in the same color, as determined by the lbl entry
%           of the current colormap.
%           If absent, then lbl = 1 is assumed.      
%
%   targetLbl only points with labels in the vector targetLbl are to
%             be plotted.  Points with other labels ignored.
%             If targetLbl = [] then all points are plotted.
%
%   plotBox  [xMin, yMin, xMax, yMax] provides the axes for
%            the plot.
%            Use plotBox = [] for automatic axis selection.
%
% Output: 
%   plotBox gives the x,y ranges of the plot (a 4-vector in the
%           same format as used by axis)
% 

if (nargin < 2) % By default, use a new figure window.
  h = figure;   
end;
if (nargin < 3) % By default, label all points with the same label.
  lbl = 1; 
end;
if (length(lbl) == 1)
  lbl = lbl * ones(size(p,1), 1);
end
if (nargin < 4) % By default, use no target labels.
  targetLbl = [];
end;
if (nargin < 5) % By default, use automatic axis ranges.
  plotBox = [];
end;

% Get the current colormap, and wrap the labels to the same range.  
cmap = colormap;
clbl = mod(lbl-1, size(colormap,1)) + 1;

% Plot each point in the corresponding color, if desired.
n = size(p,1); 
hold on;
for i=1:n
  if size(targetLbl,1) > 0    % Check for target label, plot if requested.
    if any(targetLbl==lbl(i))
      plot(p(i,1), p(i,2), 'o', 'Color', cmap(clbl(i),: ));
    end
  else                        % No specified targets, plot all points.
    plot(p(i,1), p(i,2), 'o', 'Color', cmap(clbl(i),: ));
  end
end

if size(plotBox,1) == 0       % Automatic axis adjustment.
  minPt = min(p, [], 1);
  maxPt = max(p, [], 1);
  dataRange = max(maxPt - minPt);
  minPt = minPt - 0.05 * dataRange;
  maxPt = maxPt + 0.05 * dataRange;
  plotBox = [minPt(1) maxPt(1) minPt(2) maxPt(2)];
end

% Plot with equal x,y steps, crop to the plotBox.
axis equal; axis ij;
axis(plotBox);
set(get(h,'CurrentAxes'),'Ydir','reverse');

xlabel('x'); ylabel('y');
title('Image Edgel Locations');
hold off;
