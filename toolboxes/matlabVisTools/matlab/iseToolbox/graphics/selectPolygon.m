function [selectedPoints,cornerPts] = selectPolygon(imSize,pFigure,pColor)
%
%  [selectedPoints cornerPts] = selectPolygon(imSize,pFigure,pColor)
%
%AUTHOR: Engel, Wandell
%DATE:  03.24.95
%PURPOSE:
%   This routines permits the user to select a convex polygon by
% clicking on an image.  The routine returns a binary matrix where
% the selected region is filled with 1s and the rest with 0s.  
%     1.  The region is picked according to the right hand rule (i.e.,
%         points to the right of the line are selected)
%     2.  Click left most button to select additional points, right
%         most button to end selection and close polygon
%
%ARGUMENTS:
%  imSize:  The row and column size of the image you will be selecting from
% Optional
%   pFigure: The figure number (Default is current figure)
%   pColor:  The region is outlined in this color.  Default is red.
%
%RETURNS:
% selectedPoints:  binary image of selected points
% cornerPts:       the corner points of the polygon
%

%

if nargin == 3
   pString1 = [pColor,'o'];
   pString2 = [pColor,'-'];
   figure(pFigure);
else
   pColor = 'r';
   pString1 = 'ro';
   pString2 = 'r-';
end

if nargin == 2
   figure(pFigure);
end

pts = [];
hold on
while(1)
  [x y but] = ginput(1);
  if but == 1
    pts = [pts ; x y];
  elseif but == 3
    break   
  end
  plot(pts(:,1),pts(:,2),pString1, ...
       pts(:,1),pts(:,2),pString2);
end

%	Add last line segment
nPts = size(pts,1);
line([pts(nPts,1),pts(1,1)],[pts(nPts,2),pts(1,2)],'color',pColor);

hold off

cornerPts = pts;
selectedPoints = findInteriorPoints(pts,imSize);

