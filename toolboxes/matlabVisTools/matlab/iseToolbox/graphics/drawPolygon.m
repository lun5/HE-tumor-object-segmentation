function drawPolygon(pts,pColor)
%drawPolygon(pts,pColor)
%
%AUTHOR:  Wandell
%DATE:    November 1995
%PURPOSE:
%  Early draft ... put up a polygon on the current image
%

if nargin < 2
   pColor = 'r';
end

pString1 = [pColor,'o'];
pString2 = [pColor,'-'];

nPts = size(pts,1);

for i=1:(nPts-1)
  hold on
  plot(pts(:,1),pts(:,2),pString1, ...
       pts(:,1),pts(:,2),pString2);
end

hold on
line([pts(nPts,1),pts(1,1)],[pts(nPts,2),pts(1,2)],'color',pColor);
hold off
