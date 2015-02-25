%function [fracIm, segIm, s0] = rbfFracImageNew(seed0, bias, fracDim, ...
%                                               sizeIm, useAlpha)
function [fracIm, segIm] = rbfFracImageNew(seed0, bias, fracDim, ...
                                               sizeIm, useAlpha)                                          
%function [fracIm, segIm, s0] = rbfFracImageNew(seed0, bias, fracDim, ...
%                                               sizeIm, useAlpha)
%
% seed0   is the random number seed.  (Compute new seed0 if input is <= 0)
%         (Default: 0)
% bias    The mean separation between brightnesses in the on and off
%         regions.  biases in [-4, 4] produce reasonable separations.
%         (Default: 2)
% fracDim is the fractal dimension of the texture, useable in [0 2]
%         The std dev of the grey levels in the fractal texture is 1, 
%         mean is zero.
%         (Default: 0.8)
% sizeIm  Size of output image.  The rbf's all have a fixed sigma of 15,
%         so this sets the smoothness of the boundaries.
%         (Default: [100 100]);
% useAlpha  ~=0 => use alpha map for anti-aliasing.
%
% Examples:
%   [im segIm s0] = rbfFracImage;  %% Generates new random 100x100 image
%   figure(1);clf; showIm(im); figure(2); clf; showIm(double(segIm));
%  Generates same regions and textures, but with a different bias
%   [im segIm s0] = rbfFracImage(s0, 1, 0.8);
% See also Test code at the bottom of the file

if nargin < 1 || isempty(seed0)
  seed0 = round(sum(1000*clock));
end
if nargin < 2 || isempty(bias)
  bias = 2;
end
if nargin < 3 || isempty(fracDim)
  fracDim = 0.8;
end
if nargin < 4 || isempty(sizeIm)
  sizeIm= [100 100];
end
if nargin < 5 || isempty(useAlpha)
  useAlpha = 1;
end
  
  
if seed0 <=0
  seed0 = round(sum(1000*clock));
end
%% Reset random number generator.
rand('state', seed0); randn('state', seed0);

if exist('FALSE', 'var') & ~isglobal(FALSE)
  clear FALSE TRUE
end
if ~exist('FALSE', 'var')
  global FALSE TRUE
  FALSE = (0 == 1); 
  TRUE = ~FALSE;
end

%% Set the RBF parameters
sigmaRBF = 15; sep = 1.0;  % sep is measured in sigmaRBF

%% Build Gaussian filter masks generating alpha-mask.
sigAlpha = 1.0;
sigmaSqr = sigAlpha*sigAlpha;
gFiltSize = 2 * round(3.0 * sigAlpha) + 1;
x = [1:gFiltSize] - round((gFiltSize+1)/2);
gFilt = exp(- x .* x / (2.0*sigmaSqr));
gFilt = gFilt/ sum(gFilt(:));

%%%%%%%%%%%%%%%% Generate RBF Kernels %%%%%%%%%%%%
s2 = 2*sigmaRBF^2;
nCell = round(sizeIm/(sep*sigmaRBF));
nCell = max(nCell,1);

[x y] = meshgrid(1:nCell(2), 1:nCell(1));
x = (x - 1)*(sep*sigmaRBF);
x = x + (sizeIm(2) - max(x(:)))/2;
y = (y-1) * (sep*sigmaRBF);
y = y + (sizeIm(1) - max(y(:)))/2;

ctr = [x(:) y(:)];

[x y] = meshgrid(1:sizeIm(2), 1:sizeIm(1));

%% COmpute initial gaussian kernels
im = zeros(prod(sizeIm), prod(nCell));
for k = 1:prod(nCell)
  im(:,k) = exp(-sum(([x(:) y(:)]-repmat(ctr(k,:),prod(sizeIm),1)).^2,2)/ ...
                s2);
  
end

%% COmpute normalized kernels
kern = im ./ repmat(sum(im,2),1,prod(nCell));


%%%%%%%%%%%%%%%% Generate RBF Kernels %%%%%%%%%%%%

%% Show random RBF function
c = randn(prod(nCell),1);
im = reshape(kern * c, sizeIm);

[n g] = hist(im(:), 101); 
n = n/sum(n);
c = cumsum(n);
c(end) = 1.0;

% Sample from image brightness
minArea = 0.2;
scl = (1 - minArea - minArea);
u = rand(1,1);
u = minArea + scl * u;
k = find(c>=u); k = k(1);
segIm = double(im>=g(k));

%% Compute alpha mask
if useAlpha
  alphaIm =  rconv2sep(segIm, gFilt, gFilt);
  alphaIm = min(alphaIm, 1);
  alphaIm = max(alphaIm, 0);
else
  alphaIm = segIm;
end

noise0=mkFract(sizeIm,fracDim); 
noise0 = noise0-mean(noise0(:));
noise0 = noise0 /std(noise0(:));

noise1=mkFract(sizeIm,fracDim); 
noise1 = noise1-mean(noise1(:));
noise1 = noise1 /std(noise1(:));

if FALSE
  imStats(noise0);
  imStats(noise1);
end
  
fracIm = (noise0-bias/2) .* (1 - alphaIm) + (noise1+bias/2).* alphaIm;

return;

%% Test:
s0 = round(sum(1000*clock));
for bias = -4:4
  [im segIm] = rbfFracImage(s0, bias);
  figure(1);clf; 
  subplot(1,2,1); showIm(im); 
  title(sprintf('Bias %4.1f', bias));
  subplot(1,2,2); showIm(double(segIm));
  fprintf(2, 'Press any key to continue...\n');
  pause;
end
% select bias to be 2.5 say, regenerate image
bias = 2.5;
[im segIm] = rbfFracImage(s0, bias);

figure(1);clf; 
subplot(1,2,1); showIm(im); 
title(sprintf('Bias %4.1f', bias));
subplot(1,2,2); showIm(double(segIm));
