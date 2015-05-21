function result = homogWarp(im, H, szIm, grayDefault)
 %% result = homogWarp(im, H, szIm, grayDefault)  
 %% Warp the input image by the 2D homography H.  The output
 %% image, result, is to be of size szIm, with x going from 0 to szIm(2)
 %% and y going from 0 to szIm(1).  By default, szIm = size(im).
 %% grayDefault is the gray level to use for points in the new image
 %% that cannot be mapped from the original image im.
 %% Note, the warp H is arranged such that a pixel (x,y) in 
 %% the result image comes from pixel (p1/p3, p2/p3) in the original
 %% image im, where (p1 p2 p3)^T = H * (x,y,1)^T.  (If you don't get
 %% the desired result, try using the inverse of H.)
   
 if nargin < 3
   % size of the image
   szIm = size(im);
 end
 if nargin < 4
   % default gray value
   grayDefault = 127;
 end
 
 % pixel coordinates in original image
 [x,y] = meshgrid(1:size(im,2),1:size(im,1));

 % pixel coordinates in result image
 if any(szIm ~= size(im))
   [xp,yp] = meshgrid(1:szIm(2),1:szIm(1));
   pix   = [xp(:)'; yp(:)'];
 else
   pix   = [x(:)'; y(:)'];
 end
  
% homogeneous pixels in result frame.
hPixels = [ pix; ones(1,prod(szIm))];

% corresponding warped points in original frame
hScene  = H*hPixels;

% pixel coords in original frame.
xprime=(hScene(1,:)./(hScene(3,:)))';
yprime=(hScene(2,:)./(hScene(3,:)))';
xprime = reshape(xprime, szIm);
yprime = reshape(yprime, szIm);

% Warping an image according to the homography
result = interp2(x,y,im,xprime,yprime, '*linear');
result(isnan(result)) = grayDefault;
result = reshape(result,szIm);

