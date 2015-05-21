function [rout,g,b] = imrotate(varargin)
%IMROTATE Rotate image.
%   B = IMROTATE(A,ANGLE,'method') rotates the image A by ANGLE
%   degrees.  The image returned B will, in general, be larger 
%   than A.  Invalid values on the periphery are set to one
%   for indexed images or zero for all other image types. Possible 
%   interpolation methods are 'nearest' (nearest neighbor),
%   'bilinear' or 'bicubic'.   'nearest' is the default for all 
%   images.   For indexed images, only nearest-neighbor interpolation 
%   should be used.   A can also be an M x N x 3 RGB truecolor image.
%
%   B = IMROTATE(A,ANGLE,'crop') or IMROTATE(A,ANGLE,'method','crop')
%   crops B to be the same size as A.
%
%   Without output arguments, IMROTATE(...) displays the rotated
%   image in the current axis.  
%
%   See also IMRESIZE, IMCROP, ROT90.

%   Clay M. Thompson 8-4-92
%   Copyright (c) 1992 by The MathWorks, Inc.
%   $Revision: 5.4 $  $Date: 1996/10/23 19:34:48 $

[Image,Angle,Method,ClassIn,DoCrop] = parse_inputs(varargin{:});

threeD = (ndims(Image)==3); % Determine if input includes a 3-D array

if threeD,
   r = rotateImage(Image(:,:,1),Angle,Method,DoCrop);
   g = rotateImage(Image(:,:,2),Angle,Method,DoCrop);
   b = rotateImage(Image(:,:,3),Angle,Method,DoCrop);
   if nargout==0, 
      imshow(r,g,b);
      return;
   elseif nargout==1,
      if strcmp(ClassIn,'uint8');
         rout = repmat(uint8(0),[size(r),3]);
         rout(:,:,1) = uint8(round(r*255));
         rout(:,:,2) = uint8(round(g*255));
         rout(:,:,3) = uint8(round(b*255));
      else
         rout = zeros([size(r),3]);
         rout(:,:,1) = r;
         rout(:,:,2) = g;
         rout(:,:,3) = b;
      end
   else % nargout==3
      if strcmp(ClassIn,'uint8')
         rout = uint8(round(r*255)); 
         g = uint8(round(g*255)); 
         b = uint8(round(b*255)); 
      else
         rout = r;        % g,b are already defined correctly above
      end
   end
else 
   r = rotateImage(Image,Angle,Method,DoCrop);
   if nargout==0,
      imshow(r);
      return;
   end
   if strcmp(ClassIn,'uint8')
      r = uint8(round(r*255)); 
   end
   rout = r;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: rotateImage
%

function b = rotateImage(A,ang,method,crop)
% Inputs:
%         A       Input Image
%         ang     Angle to rotate image by
%         method  'nearest','bilinear', or 'bicubic'
%         crop    1 to crop image, 0 will output entire rotated image

% Catch and speed up 90 degree rotations
if rem(ang,90)==0 & nargin<4,
   phi = round(ang - 360*floor(ang/360));
   if phi==90,
      b = rot90(A);
   elseif phi==180,
      b = rot90(A,2);
   elseif phi==270,
      b = rot90(A,-1);
   else
      b = A;
   end
   if nargout==0, imshow(b), else bout = b; end
   return
end

phi = ang*pi/180; % Convert to radians

% Rotation matrix
T = [cos(phi) -sin(phi); sin(phi) cos(phi)];

% Coordinates from center of A
[m,n,o] = size(A);
if ~crop, % Determine limits for rotated image
   siz = ceil(max(abs([(n-1)/2 -(m-1)/2;(n-1)/2 (m-1)/2]*T))/2)*2;
   uu = -siz(1):siz(1); vv = -siz(2):siz(2);
else % Cropped image
   uu = (1:n)-(n+1)/2; vv = (1:m)-(m+1)/2;
end
nu = length(uu); nv = length(vv);

blk = bestblk([nv nu]);
nblks = floor([nv nu]./blk); nrem = [nv nu] - nblks.*blk;
mblocks = nblks(1); nblocks = nblks(2);
mb = blk(1); nb = blk(2);

rows = 1:blk(1); b = zeros(nv,nu);
for i=0:mblocks,
   if i==mblocks, rows = (1:nrem(1)); end
   for j=0:nblocks,
      if j==0, cols = 1:blk(2); elseif j==nblocks, cols=(1:nrem(2)); end
      if ~isempty(rows) & ~isempty(cols)
         [u,v] = meshgrid(uu(j*nb+cols),vv(i*mb+rows));
         % Rotate points
         uv = [u(:) v(:)]*T'; % Rotate points
         u(:) = uv(:,1)+(n+1)/2; v(:) = uv(:,2)+(m+1)/2;
         if method(1)=='n', % Nearest neighbor interpolation
            b(i*mb+rows,j*nb+cols) = interp2(A,u,v,'*nearest');
         elseif all(method=='bil'), % Bilinear interpolation
            b(i*mb+rows,j*nb+cols) = interp2(A,u,v,'*linear');
         elseif all(method=='bic'), % Bicubic interpolation
            b(i*mb+rows,j*nb+cols) = interp2(A,u,v,'*cubic');
         else
            error(['Unknown interpolation method: ',method]);
         end
      end
   end
end

d = find(isnan(b));
if length(d)>0, 
   % The next line is a hack so we won't have to do the whole isind check if
   % we can see from the first 10 pixels that it is not indexed.
   if (isind(A(1:min(10,end))) & isind(A)),
      b(d) = 1; 
   else 
      b(d) = 0; 
   end
end

if isgray(A)   % This should always be true
   b = max(0,min(b,1));  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: parse_inputs
%

function [A,Angle,Method,Class,CropIt] = parse_inputs(varargin)
% Outputs:  A       the input image
%           Angle   the angle by which to rotate the input image
%           Method  interpolation method (nearest,bilinear,bicubic)
%           Class   storage class of A

A = varargin{1};
Angle = varargin{2};
Class = class(A);

switch nargin
case 2,             % imrotate        
   Method = 'nea'; 
   CropIt = 0;
case 3,
   if isstr(varargin{3}),
      Method = [lower(varargin{3}),'   ']; % Protect against short method
      Method = Method(1:3);
      if Method(1)=='c', % Crop string
         Method = 'nea';
         CropIt = 1;
      else
         CropIt = 0;
      end
   else
      error('''METHOD'' must be a string of at least three characters.');
   end
case 4,
   if isstr(varargin{3}),
      Method = [lower(varargin{3}),'   ']; % Protect against short method
      Method = Method(1:3);
   else
      error('''METHOD'' must be a string of at least three characters.');
   end
   CropIt = 1;
otherwise,
   error('Invalid input arguments');
end

if isa(A, 'uint8'),     % Convert A to Double grayscale for filtering & interpolation
   A = double(A)/255;
end

