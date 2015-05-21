%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File: oTensorDemo.m  
%%% Matlab script file. Orientation Tensor Demo
%%%
%%% Compute Canny edgels using Gaussian derivative filter kernels,
%%% and non-max suppression.  Compute orientation tensor:
%%%     oT = G(x,y; sigmaOT) *  [t(x,y) t^T(x,y)]
%%% where t(x,y) is the tangent direction for a Canny edgel (if any)
%%% at x,y.  If there is no such edgel, then t(x,y) = [0 0]^T.
%%%
%%% Dependencies
%%%   Directory: matlabVisTools
%%%      iseToolbox/pyrTools/
%%%          showIm.m  histo.m
%%%      utvisToolbox/file/
%%%          pgmRead.m resizeImageFig.m
%%%      utvisToolbox/edge/
%%%           cannyEdgels.m computeOrientTensor.m

% Things still to do:
% 1. Compare to oT = G * (gradI' * gradI)
% 2. Local max of randomness for "corner" points.

%%% ADJ, Oct. '03.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;

global matlabVisRoot

% We need to ensure the path is set for the iseToolbox.
if isempty(matlabVisRoot)
  dir = pwd;
  cd /h/51/jepson/pub/matlab   % CHANGE THIS
  startup;
  cd(dir);
end

% Check each required directories is on path (and some of the files)
which pgmRead

which cannyEdgels

FALSE = (0==1);
TRUE = ~FALSE;
SUPERIMPOSE = TRUE;
SAVE = FALSE;

% For display, scale small images up to be about 500 pixels
% in the largest dimension.  Larger images won't be scaled down.
figSz = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Blocks-World Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = pgmRead('im9h.pgm');
sclPix = round(max(10, (figSz * 10)/max(size(im))))/10;
% Display it
figure(1);clf;
showIm(im);
fprintf(2,'Press any key to continue...\n'); pause;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Other image files  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FALSE  % Select one of the following images (comment all but the
          % image you want, or copy-paste).
  
  % Noisy disk
  stepHeight = 32;
  noiseStdDev = 8.0;
  disp(['stepHeight: ' num2str(stepHeight) ...
        ', noiseStdDev: ' num2str(noiseStdDev)]);
  im = 128 + (mkDisc(64)-0.5) * stepHeight + randn(64, 64) * noiseStdDev;

  % Noisy circle
  stepHeight = 32;
  noiseStdDev = 8.0;
  im = pgmRead('circ64.pgm');
  im = 128 + (im-128)/128 * stepHeight + randn(64, 64) * noiseStdDev;

  im = pgmRead('three.pgm');
  im = pgmRead('microserf.pgm');
  im = pgmRead('im9h.pgm');
  im = pgmRead('wall.pgm');
  im = pgmRead('boxes01-add.pgm');
  im = pgmRead('parkbench.pgm');
  im = pgmRead('star_porch.pgm');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Compute Canny edgels %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 1;
minStrength = 4.0 / sigma;
[edgelIm nrmGradIm dirIm] = cannyEdgels(im, sigma, minStrength);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Draw Canny edgels %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig = figure(2); clf;
sclPix = round(max(10, 5000/max(size(im))))/10;
colormap([0 0 0; hsv(64)]);
qw = dirIm;
qw(isnan(qw)) = -0.1;
image( (ceil(qw*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);
clear qw;

% Crop whole image by default
cropBox = [1 1 size(im,2) size(im,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Crop Canny edgels %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imCrop = im(cropBox(2):cropBox(4), cropBox(1):cropBox(3));

if all(size(edgelIm) == size(im)) 
  edgelImCrop = edgelIm(cropBox(2):cropBox(4), cropBox(1):cropBox(3));
  nrmGradImCrop = nrmGradIm(cropBox(2):cropBox(4), cropBox(1):cropBox(3));
  dirImCrop = dirIm(cropBox(2):cropBox(4), cropBox(1):cropBox(3));
else
  sigma = 1;
  minStrength = 4.0 / sigma;
 [edgelImCrop nrmGradImCrop dirImCrop] = ...
     cannyEdgels(imCrop, sigma, minStrength);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Compute Orientation Tensor %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The orientation tensor is used to determine regions in which the
% edgels are oriented close to one orientation (mod pi), from other
% regions in which the edgels are more randomly oriented.
sigmaOT=3*sigma; % spatial std dev of Gaussian mask for orientation tensor
dwnSample = round(sigma); % Spatial sampling rate for orientation tensor

% Compute the orientation tensor given edgel input
orientT = computeOrientTensor(edgelImCrop, dirImCrop, sigmaOT, dwnSample);

% The orientation tensor (OT) is a 2x2 symmetric matrix at every sample point.
% The cropped edgel image has be downsampled by 
dwnSample

% So we have an image of the following size for the OT:
floor((size(edgelImCrop)-1)/dwnSample + 1)

% Since the OT is a symmetric 2x2 matrix, we need only store
% the (1,1) (1,2) and (2,2) elements.  These are stored
% in orientT(y,x, k) for k = 1, 2, 3 respectively.
size(orientT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Display Orientation Tensor %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The trace of the OT provides an image indicating the density of edgels
% in the neighbourhood
traceIm = sum(orientT(:,:,[1 3]),3);
figure(1); clf; showIm(traceIm); 
title('Trace of Orientation Tensor');
fprintf(2, 'Press any key to continue...\n');
pause;
 
% Compute the determinant of the OT (it must be >=0, except for rounding errors)
detIm = orientT(:,:,1) .* orientT(:,:,3) - orientT(:,:,2).^2;
detIm(detIm<0.0) = 0.0;

% Compute the eigenvalues of the OT (they must be >=0, except for rounding errors)
sep = traceIm .* traceIm - 4.0 * detIm;
sep(sep<0.0) = 0.0;
sep = sep.^0.5;
lam1 = (traceIm + sep)/2.0;
lam2 = (traceIm - sep)/2.0;

% The normalized difference between the eigenvalues indicates
% how uniformly directed an image region is.
directedNess = (lam1-lam2)./(traceIm + .1/(sigmaOT*sigmaOT));

% Display the directedness
figure(1); clf; showIm(directedNess .*(directedNess > 0.3));
title('Directedness of Orientation Tensor');
fprintf(2, 'Press any key to continue...\n');
pause;
 
% The normalized average of the eigenvalues indicates
% how uniformly scattered the orientations in a region are.
% This is similar to a 'corner' detector.
randomNess = 2 * lam2./(traceIm + .1/(sigmaOT*sigmaOT));
figure(1); clf; showIm(randomNess .*(randomNess > 0.3));
title('Randomness of Orientation Tensor');
fprintf(2, 'Press any key to continue...\n');
pause;
 
% Colour plot of directedness (red) and randomness (green)
% overlayed onto original image.  Notice the green pixels
% are in image regions with multiple orientations, the red
% pixels have a dominant orientation in the region.
scl = round(log(dwnSample)/log(2.0));
if scl>0
  res = blurDn(imCrop, scl);
  res = res / (2.0^scl);
else
  res = imCrop;
end
if all(size(res) == size(randomNess))
     rgb = cat(3, res*0.5 + directedNess.*(directedNess>0.3)*128, ...
                  res*0.5 + randomNess.*(randomNess>0.3)*128, ...
                  res*0.5)/255;
     rgb(rgb > 1) = 1.0;
     figure(1); clf; image(rgb);
     sclFig = max(round(figSz/max(size(randomNess))),1);
     resizeImageFig(1, size(rgb), sclFig);
     title('Directedness (r), Randomness (g)');
     fprintf(2, 'Press any key to continue...\n');
     pause;
end

%%%%%%%% Variations
% Try other images.  See commented block: Other image files 
% Select one of those images, and run remainder of this file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% END Orientation Tensor Demo %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
