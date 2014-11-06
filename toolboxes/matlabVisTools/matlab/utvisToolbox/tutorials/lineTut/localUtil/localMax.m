function [maxIm] = localMax(im, rad, val)
%% [maxIm] = localMax(im, rad, val)
%%
%% Compute the points which are local maxima of the 2D image im,
%% within a 2*rad+1 by 2*rad+1 box centered on each point.
%% The maxIm returned is equal to val at every pixel which isn't
%% a local max, and is equal to im otherwise.  
%% Optional:
%%   rad: default is rad=1.
%%   val: default is val=0.
  
if nargin < 2
  rad = 1
end
if nargin < 3
  val = 0
end

%% Protect against brain dead inputs...
if rad < 0
  rad = -rad;
elseif rad == 0
  rad = 1;
end

% create a version of the image with rad pixel wide border to make
% non-max suppression easier
[x y] = size(im);
maxIm = zeros(x+2*rad, y+2*rad);

% copy body of image to middle of 
maxIm((rad+1):(x+rad), (rad+1):(y+rad)) = im;

% duplicate border values
maxIm(1:rad, (rad+1):(y+rad)) = repmat(im(1,:), rad,1);
maxIm((x+rad+1):(x+2*rad),(rad+1):(y+rad)) = repmat(im(x,:), rad,1);
maxIm((rad+1):(x+rad),1:rad) = repmat(im(:,1), 1, rad);
maxIm((rad+1):(x+rad),(y+rad+1):(y+2*rad)) = repmat(im(:,y), 1, rad);

% Duplicate value at corners
maxIm(1:rad,1:rad) = im(1,1);
maxIm(x+rad+(1:rad),1:rad) = im(x,1);
maxIm((1:rad), y+rad+(1:rad)) = im(1,y);
maxIm(x+rad+(1:rad),y+rad+(1:rad)) = im(x,y);

% keep the value if it is at least as large as all its nbrs
newIm = maxIm;
for i=(rad+1):(x+rad)
  for j=(rad+1):(y+rad)
    if (newIm(i,j) < max(max(maxIm((i-rad):(i+rad), (j-rad):(j+rad)))))  
     newIm(i,j) = val;
    end
  end
end

% keep the important part of the image data
maxIm = newIm((rad+1):(x+rad), (rad+1):(y+rad));
