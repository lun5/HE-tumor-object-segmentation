function varargout = improfile(varargin)
%IMPROFILE Extract data along a path in an image.
%   C = IMPROFILE prompts you to define the desired path using the
%   mouse in the current figure.  The image values C along the path
%   are returned. Click at multiple points to define a multi-segment
%   path.  Left mouse button adds points to a multiple line path. 
%   Use right mouse button (shift-click on the Macintosh) or RETURN to
%   add the final point. IMPROFILE without image data only works
%   when the current axis contains an image. C = IMPROFILE(N)
%   returns N points along the selected path.  When the current axis 
%   contains an RGB truecolor image, the returned value has the R, G, 
%   and B values concatenated along the third dimension.
%
%   C = IMPROFILE(I,Xi,Yi,N) returns N image values along the 
%   path (Xi,Yi).  If N is omitted, a value for each pixel along
%   path is returned.
%
%   C = IMPROFILE(x,y,I,Xi,Yi) or C = IMPROFILE(x,y,I,Xi,Yi,N) 
%   accounts for non-default axis limits.  
%
%   [X,Y,C] = IMPROFILE(...) returns the (X,Y) position along
%   the path. PLOT3(X,Y,C) can be used to plot the results.
%
%   [X,Y,C,Xi,Yi] = IMPROFILE(...) returns the selected points in
%   (Xi,Yi). 
%
%   Without lefthand side arguments, IMPROFILE(...) plots the image
%   improfile in a separate figure.  NOTE: X and Y are returned in
%   pixel coordinates which are defined as being the center of
%   each pixel. 
%
%   Use a trailing method argument to specify the interpolation
%   method used, IMPROFILE(...,'method').  Possiblities are
%   'bilinear','bicubic', and 'nearest'.  'nearest' is the
%   default for all images.
%
%   See also IMPIXEL.

%   Copyright (c) 1993 by The MathWorks, Inc.
%   $Revision: 5.2 $  $Date: 1996/10/16 20:33:32 $

[xa,ya,a,n,method,prof,getn,getprof,u8in] = parse_inputs(varargin{:});

RGB_image = (ndims(a)==3);

% Parametric distance along segments
s = [0;cumsum(sqrt(sum((diff(prof).^2)')'))];

% Remove duplicate points if necessary.
killIdx = find(diff(s) == 0);
if (~isempty(killIdx))
   s(killIdx+1) = [];
   prof(killIdx+1,:) = [];
end

[ma,na,oa] = size(a);
xmin = min(xa(:)); ymin = min(ya(:));
xmax = max(xa(:)); ymax = max(ya(:));
dx = max( (xmax-xmin)/(na-1), eps );
dy = max( (ymax-ymin)/(ma-1), eps );

if getn,
   d = abs(diff(prof./(ones(size(prof,1),1)*[dx dy])));
   n = max(sum(max(ceil(d)')),2); % In pixel coordinates
end

% Interpolation points along segments
profi = interp1(s,prof,(0:n-1)*max(s)/(n-1));

xxa = xmin:dx:xmax;
yya = ymin:dy:ymax;

if RGB_image
   % Image values along interpolation points - r,g,b planes separately
   % Red plane
   xr = max(xxa(1),min(profi(:,1),xxa(length(xxa))));
   yr = max(yya(1),min(profi(:,2),yya(length(yya))));
   zr = interp2(xxa,yya,a(:,:,1),xr,yr,method); 
   % Green plane
   xg = max(xxa(1),min(profi(:,1),xxa(length(xxa))));
   yg = max(yya(1),min(profi(:,2),yya(length(yya))));
   zg = interp2(xxa,yya,a(:,:,2),xg,yg,method); 
   % Blue plane
   xb = max(xxa(1),min(profi(:,1),xxa(length(xxa))));
   yb = max(yya(1),min(profi(:,2),yya(length(yya))));
   zb = interp2(xxa,yya,a(:,:,3),xb,yb,method); 
else
   % Image values along interpolation points - the g stands for Grayscale
   xg = max(xxa(1),min(profi(:,1),xxa(length(xxa))));
   yg = max(yya(1),min(profi(:,2),yya(length(yya))));
   zg = interp2(xxa,yya,a,xg,yg,method);
end

if nargout == 0 % plot it
   oldgcf = gcf;
   if getprof,
      h = get(0,'children');
      fig = 0;
      for i=1:length(h),
         if strcmp(get(h(i),'name'),'Profile'),
            fig = h(i);
         end
      end
      if ~fig, % Create new window
         fig = figure('Name','Profile');
      end
      figure(fig)
   else
      fig = gcf;
   end
   if length(prof)>2
      if RGB_image
         plot3(xr,yr,zr,'r',xg,yg,zg,'g',xb,yb,zb,'b');
         set(gca,'ydir','reverse');
         xlabel X, ylabel Y;
      else
         plot3(xg,yg,zg,'b');
         set(gca,'ydir','reverse');
         xlabel X, ylabel Y;
      end
   else
      if RGB_image
         plot(sqrt((xr-xr(1)).^2+(yr-yr(1)).^2),zr,'r',...
              sqrt((xg-xg(1)).^2+(yg-yg(1)).^2),zg,'g',...
              sqrt((xb-xb(1)).^2+(yb-yb(1)).^2),zb,'b');
         xlabel('Distance along profile');
      else
         plot(sqrt((xg-xg(1)).^2+(yg-yg(1)).^2),zg,'b');
         xlabel('Distance along profile');
      end
   end
else
   x = xg;
   y = yg;
   if RGB_image
      zg = cat(3,zr,zg,zb);
   end
   if isa(zg,'double') & u8in   %Output the same data type we got in
      z = uint8(round(zg*255));
   else
      z = zg;
   end
   xi = prof(:,1);
   yi = prof(:,2);
   switch nargout
   case 1,
      varargout{1} = z;
   case 3,
      varargout{1} = x;
      varargout{2} = y;
      varargout{3} = z;
   case 5,
      varargout{1} = x;
      varargout{2} = y;
      varargout{3} = z;
      varargout{4} = xi;
      varargout{5} = yi;
   otherwise
      error('Invalid output arguments.');
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: parse_inputs
%

function [Xa,Ya,Img,N,Method,Prof,GetN,GetProf,u8in]=parse_inputs(varargin)
% Outputs:
%     Xa        2 element vector for non-standard axes limits
%     Ya        2 element vector for non-standard axes limits
%     A         Image Data
%     N         number of image values along the path (Xi,Yi) to return
%     Method    Interpolation method: 'nearest','bilinear', or 'bicubic'
%     Prof      Profile Indices
%     GetN      Determine number of points from profile if true.
%     GetProf   Get profile from user via mouse if true also get data from image.
%     u8in      True when the input is a uint8 image

% Set defaults
N = [];
GetN = 1;    
GetProf = 0; 
GetCoords = 1;  %     GetCoords - Determine axis coordinates if true.
 
Method = 'nearest';

switch nargin
case 0,            % improfile
   GetProf = 1; 
   GetCoords = 0;
   
case 1,            % improfile(n) or improfile('Method')
   if isstr(varargin{1})
      Method = varargin{1}; 
   else 
      N = varargin{1}; 
      GetN = 0; 
   end
   GetProf = 1; 
   GetCoords = 0;
   
case 2,            % improfile(n,'method')
   if isstr(varargin{2}),
      Method = varargin{2};
      N = varargin{1}; 
      GetN = 0; 
   else 
      error('Wrong number of arguments or unknown interpolation method.');
   end
   GetProf = 1; 
   GetCoords = 0;
   
case 3,   % improfile(a,xi,yi)
   A = varargin{1};
   Xi = varargin{2}; 
   Yi = varargin{3}; 
   
case 4,   % improfile(a,xi,yi,n) or improfile(a,xi,yi,'method')
   A = varargin{1};
   Xi = varargin{2}; 
   Yi = varargin{3}; 
   if isstr(varargin{4}) 
      Method = varargin{4}; 
   else 
      N = varargin{4}; 
      GetN = 0; 
   end
   
case 5, % improfile(x,y,a,xi,yi) or improfile(a,xi,yi,n,'method')
   if isstr(varargin{5}), 
      A = varargin{1};
      Xi = varargin{2}; 
      Yi = varargin{3}; 
      N = varargin{4}; 
      Method = varargin{5}; 
      GetN = 0; 
   else
      GetCoords = 0;
      Xa = varargin{1}; 
      Ya = varargin{2}; 
      A = varargin{3};
      Xi = varargin{4}; 
      Yi = varargin{5}; 
   end
   
case 6, % improfile(x,y,a,xi,yi,n) or improfile(x,y,a,xi,yi,'method')
   Xa = varargin{1}; 
   Ya = varargin{2}; 
   A = varargin{3};
   Xi = varargin{4}; 
   Yi = varargin{5}; 
   if isstr(varargin{6}), 
      Method = varargin{6}; 
   else 
      N = varargin{6};
      GetN = 0; 
   end
   GetCoords = 0;
   
case 7, % improfile(x,y,a,xi,yi,n,'method')
   if ~isstr(varargin{7}), 
      error('Method must be ''bilinear'',''bicubic'', or ''nearest''.');
   end
   Xa = varargin{1}; 
   Ya = varargin{2}; 
   A = varargin{3};
   Xi = varargin{4}; 
   Yi = varargin{5}; 
   N = varargin{6};
   Method = varargin{7}; 
   GetN = 0;
   GetCoords = 0; 
   
otherwise
   error('Invalid input arguments.');
end

if ~GetN,
   if N<2, error('N must be > 1.'); end
end

if GetCoords & ~GetProf,
   Xa = [1 size(A,2)];
   Ya = [1 size(A,1)];
end

if GetProf, % Get profile from user if necessary using data from image
   [Xa,Ya,A,state] = getimage;
   if prod(size(A))==1, 
      error('Can''t get profile of 1-by-1 image.'); 
   end
   if ~state, 
      error('Requires an image in the current axis.'); 
   end
   Prof = getline(gcf); % Get profile from user
else  % We already have A, Xi, and Yi
   if prod(size(A))==1, 
      error('Can''t get profile of 1-by-1 image.'); 
   end
   if prod(size(Xi))~=prod(size(Yi)),
      error('Xi and Yi must have the same number of points');
   end
   Prof = [Xi(:) Yi(:)]; % [xi yi]
end

u8in = isa(A,'uint8');

% Promote the image to double if we aren't using nearest
if u8in & ~strcmp(Method,'nearest')
   Img = double(A)/255;
else
   Img = A;
end

   