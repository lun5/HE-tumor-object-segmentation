function [selectedPts, endPts, lHandle]= selectLine(isize,thickness,overlay)
% 
% [selectedPts, endPts, lHandle] = selectLine(isize,thickness,overlay)
%
%AUTHOR:  Wandell
%DATE:    April 26, 1995
%PURPOSE:  
%
%  Select a line from image data displayed in the axis('image') format.
%
%	isize:  The size of the image that line is selected from
%	thickness:  The thickness of the selected region (default = 1)
%	overlay:    Draw an overlay line (1 = true, 0 = default)
%
%	selectedPts:  An array if isize with 1s in the selected positions
%		and zeros elsewhere.
% 	endPts: The locations of the line endpoints.  To draw an overlay 
%		line use:
%	lHandle: line(endPts(:,1),endPts(:,2),'color','k');
%

if nargin < 3
  overlay = 0;
end
if nargin < 2
  thickness = 1;
end

% First column is horizontal scrren, 2nd vertical
endPts = fix(ginput(2));

x1 = endPts(1,1); y1 = endPts(1,2);
x2 = endPts(2,1); y2 = endPts(2,2);
y = [y1; y2];
x = [x1, 1; x2 , 1];
selectedPts = zeros(isize);

if overlay == 1
 X = [x1 ; x2];
 Y = [y1 ; y2];
 lHandle = line(X,Y,'color','r');
end

if thickness == 1
 disp('Thin selection')
 if x1 ~= x2	
%
% Compute the parameters of the line y = ax + b
%
  ab = inv(x)*y;
  slope = ab(1);
  offset = ab(2);

%
% Fill in the values along the longest dimension
%
  if abs(x2 - x1) > abs(y2 - y1)
   xValues = min(x1,x2):max(x1,x2);
   yValues = slope*xValues + offset;
  else
   yValues = min(y1,y2):max(y1,y2);
   xValues = (yValues - offset)/slope;
  end

 else  %The equation is x = constant;
  yValues = min(y1,y2):max(y1,y2);
  xValues = ones(size(yValues))*x1;
 
 end
 selectedPts =  ...
   setMatrixEntries(xValues,yValues, selectedPts,ones(size(yValues)));

else

 disp('Thickness selection')
 pt1 = [x1 , y1]; pt2 = [x2 , y2];
 selectedPts = perpDistance(pt1,pt2,isize);
 selectedPts = reshape(selectedPts,isize(1),isize(2));

 y = min(endPts(:,1)):max(endPts(:,1));
 x = min(endPts(:,2)):max(endPts(:,2));
 box = zeros(size(selectedPts));
 box(x,y) = ones(size(box(x,y)));

%
%  Each position on the line must less than 1/sqrt(2) pixel
%  away from precisely one, except if it perfectly splits the box.
%  Occasional bad luck, hunh.  
%  This does select jagged lines.
 selectedPts = (abs(selectedPts) < 0.71*thickness) & (box == 1);

end

