function dist = perpDistance(pt1,pt2,imSize);
%
%   dist = perpDistance(pt1,pt2,imSize);
%
%AUTHOR:  Engel (mostly), Wandell
%DATE:    April, 1995
%PURPOSE:
%
% Takes two points and returns a vector of distances (with left-hand
% rule sign) of the distances to the line from each point in the image
%	pt1, pt2:  (x,y) coordinates of the two points
%	imSize:    image size
%
%	dist:	   the distance of each of the points in
%		   an array of imSize from the line defined by the points
%

%pt1 = [0 0];
%pt2 = [1 1];
%imSize = [10 10];

a = pt1(1); b = pt1(2);
dx = pt2(1)-pt1(1);
dy = pt2(2)-pt1(2);
d = [dx dy]';
d = d/sqrt(dx*dx + dy*dy); 
[u s v] = svd(d');
perp = v(:,2);
perp = makeLeftHandVector(d,perp);

%
%  Create all the vector locations
[xcoord,ycoord] = meshgrid(1:imSize(2),1:imSize(1));

%  Shift the image points down so that a,b is the origin
%
u = a-xcoord(:);                         % Vector from start of line 
v = b-ycoord(:);

%
% Compute distance for each point in image to the line perpendicular
% to the selected line. 
%


dist = [u , v]*perp;
