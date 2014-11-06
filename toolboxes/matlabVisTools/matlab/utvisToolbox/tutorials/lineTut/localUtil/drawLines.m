function [] = drawLines(par, cropBox)
% drawLines(par, cropBox)
% Draw lines in the current figure given the line
% parameters for the k-th line are:
%     par(1,k) x + par(2,k) y + par(3,k) = 0

  hold on;
  for j = 1:size(par,1)
    endPts = cropLineInBox(par(j,1:2), par(j,3), cropBox);
    if ~isnan(endPts(1,1))
      % Covering line
      plot(endPts(:,1), endPts(:,2),'r-', 'LineWidth', 0.5);
    end
  end
  hold off;