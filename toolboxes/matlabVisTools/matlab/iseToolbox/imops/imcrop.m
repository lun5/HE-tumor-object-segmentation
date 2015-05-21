function [bout,rect] = imcrop(x,y,a,rect)
%IMCROP Crop image.
%   B = IMCROP prompts you to define a rectangular subimage
%   of the current figure using the mouse and returns in B the
%   interior of the bounded area.  IMCROP only works if the 
%   current axis contains an image (i.e., created using IMSHOW,
%        IMAGE, or, IMAGESC).  Example:
%           imshow(X,map) % Display image
%           b = imcrop;
%
%   B = IMCROP(A,RECT) returns the submatrix with the bounding
%   rectangle RECT = [xmin, ymin, width, height] in spatial
%   coordinates, where x increases from left to right, and y 
%   increases from top to bottom.  (i.e. [.5 .5 N M] for an 
%   M-by-N image). 
%   
%   B = IMCROP(x,y,A,RECT) accounts for non-default axis limits.
%
%   [B,RECT] = IMCROP(...) returns the selection rectangle.
%   IMCROP without output arguments displays the cropped image.
%
%   See also IMZOOM.

%   Clay M. Thompson 1-22-93
%   Copyright (c) 1993 by The MathWorks, Inc.
%   $Revision: 5.3 $  $Date: 1996/10/24 14:50:30 $

if nargin==0, % Get information from the current figure
   [x,y,a,hasimage] = getimage;
   if ~hasimage,
      error('The current figure must contain an image to use IMCROP.');
   end
   rect = getrect(gcf); % Get rect info from the user.
elseif nargin==1 | nargin==3,
   error('Wrong number of input arguments.');
elseif nargin==2,
   a = x;
   [m,n,o] = size(a);
   rect = y;
   x = 1:n; y = 1:m;
end

[m,n,o] = size(a);
xmin = min(x(:)); ymin = min(y(:));
xmax = max(x(:)); ymax = max(y(:));
% Transform rectangle into row and column indices.
kx = n-1; ky = m-1;
i1 = max(ceil((rect(2)-ymin)/(ymax-ymin)*ky+1),1);
j1 = max(ceil((rect(1)-xmin)/(xmax-xmin)*kx+1),1);
i2 = min(floor((rect(2)+rect(4)-ymin)/(ymax-ymin)*ky+1),m);
j2 = min(floor((rect(1)+rect(3)-xmin)/(xmax-xmin)*kx+1),n);

b = a(i1:i2,j1:j2,:);

if nargout==0,
   imshow(b)
   return
end
bout = b;
