function [orientT] = computeOrientTensor(edgelIm, dirIm, sigmaOT, dwnSamp)
% [orientT] = computeOrientTensor(edgelIm, dirIm, sigmaOT, dwnSamp)
%
% Compute the orientation tensor given edgels using  a Gaussian filter.
%
% PARAMS:  
%  IM the 2D gray level input image
%  SIGMAOT (optional, default = 4.0) the std dev of the Gaussian filter.
%
% OUTPUT:  
%  orientT  orientation tensor  [G conv tx*tx, G conv tx*ty, G conv ty * ty]

% ADJ 9/01.

%%%%%%%%%%%%%%% Check parameters  %%%%%%%%%%%%%%%%%%%%%%%

%%% Fill in default parameter value, and correct bogus parameter values.
if ~exist('sigmaOT', 'var') % sigmaOT is the std dev for the Gaussian filters
  sigmaOT = 1.0;
elseif (sigmaOT < 0.5/3)
  fprintf(2, 'User specified sigmaOT = %e is too small.\n', sigmaOT); 
  sigmaOT = 0.5000001/3;
  fprintf(2, 'Using sigmaOT = %f instead\n', sigmaOT);
end
if ~exist('dwnSamp', 'var') % sigmaOT is the std dev for the Gaussian filters
  dwnSamp = round(sigmaOT);
end
if (dwnSamp < 1)
  dwnSamp = 1;
end


%%%%%%%%%%%%%%% Build Filter Kernels %%%%%%%%%%%%%%%%%%%%%%%

%%% Build Gaussian filter masks.
sigmaSqr = sigmaOT*sigmaOT;
gFiltSize = 2 * round(3.0 * sigmaOT) + 1;
x = [1:gFiltSize] - round((gFiltSize+1)/2);
gFilt = exp(- x .* x / (2.0*sigmaSqr));
gFilt = gFilt/ sum(gFilt(:));

%%% Build tangent image
tmpIm = zeros(size(edgelIm));
tangIm = zeros([size(edgelIm), 2]);
tmpIm(edgelIm) = -sin(pi*dirIm(edgelIm));
tangIm(:,:,1) = tmpIm;
tmpIm = zeros(size(edgelIm));
tmpIm(edgelIm)  = cos(pi*dirIm(edgelIm));
tangIm(:,:,2) = tmpIm;
clear tmpIm;

%%%%%%%%%%%%%%% Do separable convolutions %%%%%%%%%%%%%%%%%%%
dwnSize = floor((size(edgelIm) - 1)/dwnSamp) + 1;

%% Column filtering
gFilt = gFilt(:);

orientT = zeros([dwnSize, 3]);
tmpIm = corrDn(tangIm(:,:,1) .* tangIm(:,:,1), gFilt, ...
               'reflect1', [dwnSamp 1]);
orientT(:,:,1) = corrDn(tmpIm, gFilt', 'reflect1', [1 dwnSamp]);

tmpIm = corrDn(tangIm(:,:,1) .* tangIm(:,:,2), gFilt, ...
               'reflect1', [dwnSamp 1]);
orientT(:,:,2) = corrDn(tmpIm, gFilt', 'reflect1', [1 dwnSamp]);

tmpIm = corrDn(tangIm(:,:,2) .* tangIm(:,:,2), gFilt, ...
               'reflect1', [dwnSamp 1]);
orientT(:,:,3) = corrDn(tmpIm, gFilt', 'reflect1', [1 dwnSamp]);

