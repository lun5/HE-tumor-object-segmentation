function [result,A] = imWarpAffine(im,A, crop,method)
%
% function [result,mA]= imWarpAffine(im,A, crop)
%
% im: input image
% A: 2x3 affine transform matrix, eg [R [t_x; t_y]]
% crop = 0 indicates the output result should be resized to
%          contain the warped image (using only integer shifts).
%      = 1 use original image coordinates to crop warped image
% method: interpolation method (eg 'bilinear', default = 'bicubic').
%
% if a transformed point is outside of the volume, NaN is used
% also returns the modified affine warp (from deformed to original image), mA
%
% Allan Jepson
% May 2003

if (size(A,1)>=2)
  A=[A(1:2,:); 0 0 1];
end

if ~exist('method')
   method = 'bicubic';
end

%% Compute bounding box of warped image
if crop
  wIm = zeros(size(im));
else
  corners = [1 1 size(im,2) size(im,2); 1 size(im,1) 1 size(im,1)];
  inA = inv(A);
  warpedCoords = inA * [ corners; ones(1, size(corners,2))];
  xmin = floor(min(warpedCoords(1,:))+100*eps);
  xmax = ceil(max(warpedCoords(1,:))-100*eps);
  xmax = max(xmin+1, xmax);
  ymin = floor(min(warpedCoords(2,:))+100*eps);
  ymax = ceil(max(warpedCoords(2,:))-100*eps);
  ymax = max(ymin+1, ymax);
  
  wIm = zeros(ymax-ymin+1, xmax-xmin+1);
  
  inA(1:2,3) = inA(1:2,3) + [1 - xmin; 1-ymin];
  A(1:2, 3) = -A(1:2, 1:2) * inA(1:2,3);
end
  

% Compute coordinates corresponding to input 
% and transformed coordinates for result
[x,y]=meshgrid(1:size(wIm,2),1:size(wIm,1));
coords=[x(:)'; y(:)'];
homogeneousCoords=[coords; ones(1,prod(size(wIm)))];
warpedCoords=A*homogeneousCoords;

xprime=warpedCoords(1,:)';
xprime(xprime<1.0 & xprime > 1.0-10*eps) = 1.0;
xprime(xprime>size(im,2) & xprime < size(im,2)+10*eps) = size(im,2);

yprime=warpedCoords(2,:)';
yprime(yprime<1.0 & yprime > 1.0-10*eps) = 1.0;
yprime(yprime>size(im,1) & yprime < size(im,1)+10*eps) = size(im,1);

[x,y]=meshgrid(1:size(im,2),1:size(im,1));
%result = interp2(x,y,im,reshape(xprime, size(wIm)),reshape(yprime,
%size(wIm)));
result = interp2(x,y,im, xprime, yprime,method);
result = reshape(result,size(wIm));

return;



%%% Debug

%im=[1 2 3; 4 5 6; 7 8 9; 10 11 12]';
%A= [0 1  0;
%    -1 0 0;
%    0 0 1];
%A= [1 0 0;
%    0 1 .5;
%    0 0 1];
%A= [1 0 .5;
%    0 1 .5;
%    0 0 1];
%res=imWarpAffine(im,A,0)
%im = pgmRead('einstein.pgm');
%im = im(32:(end-32),:);
%size(im)
%A= [0 1  0;
%    -1 0 0;
%    0 0 1];
%
%tmp = im;
%for k = 1:4
%  res = imWarpAffine(tmp,A,0);
%  size(res)
%  showIm(res);
%  pause;
%  tmp = res;
%end
%showIm(im-res)



