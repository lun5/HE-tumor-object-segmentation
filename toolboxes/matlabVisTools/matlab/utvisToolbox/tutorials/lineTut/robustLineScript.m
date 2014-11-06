%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File: lineTutorial.m  
%%% Matlab script file
%%%
%%% Compute Canny edgels using Gaussian derivative filter kernels,
%%% and non-max suppression.  Fit lines to edgel data.

%%% Dependencies
%%%   Directory: matlabVisTools
%%%      iseToolbox/pyrTools/
%%%           mkDisc.m showIm.m  histo.m
%%%      utvisToolbox/file/
%%%          pgmRead.m resizeImageFig.m
%%%      utvisToolbox/edge/
%%%           cannyEdgels.m computeOrientTensor.m
%%%      ./localUtil/

% Things still to do:
% 1. Robust line estimator does not have an associated tutorial, just
%    the Matlab code.
% 2. Seeds from connected groups of edgels.  Lowe's algorithm.
% 3. Quadratic and cubic curve segments.

% ADJ, Oct. '01.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;

global matlabVisRoot
% We need to ensure the path is set for the iseToolbox and utvisToolbox.
if isempty(matlabVisRoot)
  dir = pwd;
  cd ~jepson/pub/matlab   % CHANGE THIS
  startup;
  cd(dir);
end

lineRoot = [matlabVisRoot '/utvisToolbox/tutorials/lineTut']; 
addpath([lineRoot '/localUtil']);

% Check each required directories is on path (and some of the files)
which mkDisc
which pgmRead
which cannyEdgels

FALSE = (0==1);
TRUE = ~FALSE;
SUPERIMPOSE = TRUE;
SAVE = FALSE;

% Random number generator seed:
seed = round(sum(1000*clock));
seed0 = seed;  % Save the seed for restarts.
% Start the random number generators rand given this seed.
rand('state', seed);
% We also need to start randn using the autmatically generated seedn:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

% For display, scale small images up to be about 500 pixels
% in the largest dimension.  Larger images won't be scaled down.
figSz = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Blocks-World Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = pgmRead('im9h.pgm');
sclPix = round(max(10, (figSz * 10)/max(size(im))))/10;
% Display it
figure(1);clf;
showIm(im);
fprintf(2,'Press any key to continue...\n'); pause;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Robust line segment estimation  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FALSE  % Select one of the following images.
  
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
  im = pgmRead('parkbench.pgm');
  im = pgmRead('star_porch.pgm');
  im = pgmRead('clutter.pgm');
  imPath = 'c:/images/boxes/'
  fname = 'boxes01-add';
  fname = 'boxes02-add';
  fname = 'boxes03-add';
  fname = 'boxes04-add';
  fprintf(2, 'Reading %s\n', [imPath fname '.pgm']);
  im = pgmRead([imPath fname '.pgm']);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Compute Canny edgels %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 1;
minStrength = 4.0 / sigma;
[edgelIm nrmGradIm dirIm] = cannyEdgels(im, sigma, minStrength);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Draw Canny edgels %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig = figure(2); clf;
sclPix = round(max(10, 5000/max(size(im))))/10;
colormap([0 0 0; hsv(64)]);
qw = dirIm;
qw(isnan(qw)) = -0.1;
image( (ceil(qw*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);
clear qw;

% Crop whole image by default (or use the above code for mousing
% in a crop box).
cropBox = [1 1 size(im,2) size(im,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Crop Canny edgels %%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%% Compute Orientation Tensor %%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%% Display Orientation Tensor %%%%%%%%%%%%%%%%%%%%%%
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
directedNess = (lam1-lam2)./(traceIm + .3/(sigmaOT*sigmaOT));

% Display the directedness
figure(1); clf; showIm(directedNess .*(directedNess > 0.3));
title('Directedness of Orientation Tensor');
fprintf(2, 'Press any key to continue...\n');
pause;
 
% The normalized average of the eigenvalues indicates
% how uniformly scattered the orientations in a region are.
% This is similar to a 'corner' detector.
randomNess = 2 * lam2./(traceIm + .3/(sigmaOT*sigmaOT));
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
     resizeImageFig(1, size(rgb), round(2^scl));
     title('Directedness (r), Randomness (g)');
     fprintf(2, 'Press any key to continue...\n');
     pause;
end

% Clean up orientation tensor stuff (except directedNess, which
% we use below.
clear randomNess traceIm sep detIm rgb res orientT lam1 lam2
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Associate the edgels with the directedNess image %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the edgels and their orientations
[p theta] = getEdgelData(edgelImCrop, dirImCrop);

% Build indices from edgels into the directedNess array,
% accounting for the downsampling.
dwnIndex = cat(2, floor((p(:,2)-1)/dwnSample + 1), ...
                  floor((p(:,1)-1)/dwnSample + 1));
dwnIndex = dwnIndex * [1; size(directedNess, 1)] - size(directedNess,1);

% For line estimation intial guesses, only use seeds with
% directedNess larger than 0.35.
directedSeeds = (directedNess(dwnIndex) > 0.35);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Robust estimation of lines   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
lines = robustLineEstimator(p, theta, sigma, imCrop, directedSeeds);
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Display estimated lines   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

hFig = figure(2);
clf;
if SUPERIMPOSE
  image(imCrop);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imCrop), 2); hold on;
axis equal;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
% hp = plot(p(:,1), p(:,2), 'or', 'MarkerSize', 1.5, 'MarkerFaceColor', 'r');

for k=1:size(lines,2)
  % Plot k-th line from x0,y0 to x1,y1 
  plot(lines([1 3],k), lines([2 4],k),'g-', 'LineWidth', 2.0);
end

% print -depsc -r600 microserf_overlay_lines.eps
% print -depsc -r600 microserf_im.eps
% print -depsc -r600 star_porch_overlay_lines.eps
% print -depsc -r600 star_porch_im.eps
% print -depsc -r600 parkbench_overlay_lines.eps
% print -depsc -r600 parkbench_im.eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Save estimated lines   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if SAVE

  fprintf(2, 'Writing %s\n', [imPath fname '-lines' '.pfm']);
  pfmWrite(lines, [imPath fname '-lines' '.pfm']);
  
end
