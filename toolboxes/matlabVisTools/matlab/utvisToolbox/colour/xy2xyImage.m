function xyImCoords = xy2xyImage(xyCoords)
% function xyImCoords = xy2xyImage(xyCoords)
%
% Convert the xyCoords in the columns of xyCoords to
% the corresponding pixel in the picture of the xy colour
% space, diagxy.tif.

% Data from moused in points:
%%  (x,y) = (0, 0.1) at pixel  (px,py) = (51, 290)
%%  (x,y) = (0.7, 0.8) at pixel  (px,py) = (323, 21)
  dPix = [(323 - 51)/0.7; (21 - 290)/0.7 ];  
  p0 = [51; 290 - dPix(2) * 0.1];
  xyImCoords = repmat(p0, [1 size(xyCoords,2)]) + diag(dPix) * xyCoords;
