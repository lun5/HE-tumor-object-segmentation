%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File: cannyTutorial.m  
%%% Matlab script file
%%%
%%% Compute Canny edgels using Gaussian derivative filter kernels,
%%% and non-max suppression.  Try it on several images.

%%% Hysteresis thresholding and sub-pixel resolution are NOT
%%% implemented.

%%% Dependencies
%%%   Directory:
%%%      iseToolbox/pyrTools/
%%%           mkDisc.m showIm.m  histo.m
%%%      utvisToolbox/file/
%%%          pgmRead.m resizeImageFig.m
%%%      utvisToolbox/tutorials/
%%%          cannyEdgelDemo.m

%% Things still to do:
%%  Compare with Robert's and Sobel edgel operators.
%%  Sub-pixel edgel resolution.
%%  Comparison with G2/H2 edgels.

%%% ADJ, Oct. '01.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;

%% Set the root for the vision toolbox.
%% This is already done if you ran the matlabVisRoot startup file.
global matlabVisRoot
if size(matlabVisRoot,1) == 0
  matlabVisRoot = '/h/u1/jepson/pub/matlab';
end

%% Set the root for the different toolbox directories.
iseRoot = [matlabVisRoot '/iseToolbox'];
utvisRoot = [matlabVisRoot '/utvisToolbox'];
imageRoot = [matlabVisRoot '/images'];

%% Check each required directories is on path (and some of the files)
which mkDisc
which pgmRead
which cannyEdgelDemo

%% If these files cannot be found, then update your path
path(path, [iseRoot '/pyrTools/MEX']);
path(path, [iseRoot '/pyrTools']);
path(path, [utvisRoot '/file']);
path(path, [utvisRoot '/tutorials']);

%% Set the logical constants
global FALSE
global TRUE
FALSE = 0;
TRUE = ~FALSE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Noisy Disk Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate noisy disk image
stepHeight = 32;
noiseStdDev = 8.0;
disp(['stepHeight: ' num2str(stepHeight) ...
      ', noiseStdDev: ' num2str(noiseStdDev)]);
im = 128 + (mkDisc(64)-0.5) * stepHeight + randn(64, 64) * noiseStdDev;

%% Display it
figure(1);clf;
showIm(im);
title(sprintf('Noisy disk (SNR = %f)', stepHeight/noiseStdDev));
fprintf(2,'Noisy disk image (SNR = %f/%f)\n', stepHeight, noiseStdDev);
fprintf(2,'Press any key to continue...\n'); pause;

%% Run Canny edges on this image, and show the steps.
%% You should open the file cannyEdgelDemo.m in an editor
%% and follow along in the code during the execution.  The comments
%% explain what is being done.
sigma = 2;
hideSteps = FALSE;  
cannyEdgelDemo(im, sigma, 0.0, hideSteps);
 
%% There are several alternative ways to run the cannyEdgelDemo
help cannyEdgelDemo;

%% Quiet version.
hideSteps = TRUE;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma);
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
showIm(dirIm);
fprintf(2,'Press any key to continue...\n'); pause;

%% Show same edgel image with a periodic colour map
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);
fprintf(2,'Press any key to continue...\n'); pause;

%% Use prespecified edge amplitude threshold
minStrength = 3.0;
fprintf(2,'Recompute edgels with given amplitude threshold %f\n',...
        minStrength );
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, minStrength);

%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);

%% Conclude: 
%%   Edgels and their orientations seem stable given image noise.
%%   Non-max suppression is sensistive to choice of threshold.
fprintf(2,'Press any key to continue...\n'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Noisy Circle Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate noisy circle image
stepHeight = 32;
noiseStdDev = 8.0;
disp(['stepHeight: ' num2str(stepHeight) ...
      ', noiseStdDev: ' num2str(noiseStdDev)]);
im = pgmRead('circ64.pgm');
im = 128 + (im/255.0-0.5) * stepHeight + randn(64, 64) * noiseStdDev;

%% Display it
figure(1);clf;
showIm(im);
title(sprintf('Noisy circle (SNR = %f)', stepHeight/noiseStdDev));
fprintf(2,'Noisy circle image (SNR = %f/%f)\n', stepHeight, noiseStdDev);
fprintf(2,'Press any key to continue...\n'); pause;

%% Run Canny edges on this image, and show the steps.
sigma = 2;
hideSteps = FALSE;  
cannyEdgelDemo(im, sigma, 0.0, hideSteps);

%% Notice the edgels straddle the circle in the image.
 
% Run Canny edges at a finer scale on this image, and show the steps.
sigma = 1;
hideSteps = FALSE;  
cannyEdgelDemo(im, sigma, 0.0, hideSteps);
 
%% Conclude: 
%%   Edgels straddle narrow image curves.
fprintf(2,'Press any key to continue...\n'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Handwriting Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = pgmRead('three.pgm');
%% Display it
figure(1);clf;
showIm(im);
fprintf(2,'Press any key to continue...\n'); pause;

% Compute Canny edgels
hideSteps = FALSE;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, 0.0, hideSteps);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);

%% Conclude: 
%%   Gradient amplitude histogram can provide a good clue to 
%%   an appropriate choice for the gradient amplitude threshold.
fprintf(2,'Press any key to continue...\n'); pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Blocks-World Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = pgmRead('im9h.pgm');
%% Display it
figure(2);clf;
showIm(im);
fprintf(2,'Press any key to continue...\n'); pause;

% Compute Canny edgels
sigma = 1;
hideSteps = FALSE;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma);

%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);

fprintf(2,'Read the text of cannyTutorial.m for more info here\n');
%% Try to lower the gradient amplitude threshold to get
%% edgels along the center edge of the leftmost block in the image.
%% It probably has relatively few edgels on it in your first try.
%% (Simply repeat the call to cannyEdgelDemo above, and mouse
%% in a smaller threshold).

%% Also, note the double edgels on the right edge of the top cube,
%% and the left side of the rightmost cube.
%% Blow up figure(2) (the original image) in these neighbourhoods.
%% Note the extra brightness changes near these edges.  These
%% may be image compression artifacts (the original image was
%% captured and compressed using JPEG).  However, highlights
%% also often appear on or near sharp folds in the surface. And
%% these can produce a similar effect.

%% Conclude: 
%%   Common errors in detecting edgels, even in simple images, include:
%%     - false positives (extra edgels due to noise/shadows/texture).
%%     - false negatives (missing edgels due to low signal-to-noise).     
%%     - double edges can appear near sharp contrast changes.
fprintf(2,'Press any key to continue...\n'); pause;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Einstein Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = pgmRead('einstein.pgm');
%% Display it
figure(1);clf;
showIm(im);
fprintf(2,'Press any key to continue...\n'); pause;

%% Note: the amplitude histogram for this image is more difficult
%% to interpret (i.e. to set an appropriate minStrength threshold).
%% The amplitude threshold can be set pretty low.
sigma = 2.0;
hideSteps = TRUE;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, 0.0, hideSteps);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);
fprintf(2,'Press any key to continue...\n'); pause;

% Try to resolve more detail by using a smaller sigma
sigma = 1.0;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, 0.0, hideSteps);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);

%% Notice the increased detail in the hair and coat.
fprintf(2,'Press any key to continue...\n'); pause;

%% Check stability of edgels to additive (independent) image noise.
nStdDev = 10.0;
noisyIm = im + randn(size(im)) * nStdDev;
noisyIm(noisyIm < 0.0) = 0.0;
noisyIm(noisyIm > 255.0) = 255.0;

% display noisy image
hFig = figure(1); clf;
image(noisyIm);
colormap(gray(256));
resizeImageFig(hFig, size(im), sclPix);
fprintf(2,'Press any key to continue...\n'); pause;

% Run edgel detector
sigma = 2.0;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(noisyIm, sigma, 0.0, hideSteps);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);
fprintf(2,'Press any key to continue...\n'); pause;

% Try to resolve more detail by rerunning with: 
sigma = 1.0 
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(noisyIm, sigma, 0.0, hideSteps);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);

%% Conclude: 
%%   Smaller values of sigma generate edgels for finer detail,
%%   but are sensitive to image noise.
%%   Gradient amplitude histogram does not always provide a clue to 
%%   an appropriate gradient amplitude threshold.
fprintf(2,'Press any key to continue...\n'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Microserf Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

im = pgmRead('microserf.pgm');
%% Display it
figure(1);clf;
showIm(im);

fprintf(2,'Press any key to continue...\n'); pause;

hideSteps = TRUE;
sigma = 2;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, 0.0, hideSteps);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(2); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);
% Blow up the word 'noise' in the bottom right corner of this figure.

%% Try a finer resolution
sigma = 1.0;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, 0.0, hideSteps);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(3); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);
% Blow up the word 'noise' in the bottom right corner of this
% figure.  Compare with figure(2);   Note the i,s,e in 'noise'
% are cleaner in the edgel image formed with sigma = 1.

%% Check stability of edgels to additive (independent) image noise.
nStdDev = 20.0;
noisyIm = im + randn(size(im)) * nStdDev;
noisyIm(noisyIm < 0.0) = 0.0;
noisyIm(noisyIm > 255.0) = 255.0;

hFig = figure(1); clf;
image(noisyIm);
colormap(gray(256));
resizeImageFig(hFig, size(im), sclPix);
fprintf(2,'Press any key to continue...\n'); pause;

%% Get edgels of noisy image.
sigma = 2.0;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(noisyIm, sigma);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(4); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);

%% Try a finer resolution
sigma = 1.0;
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(noisyIm, sigma);
%% Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(5); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);
% Blow up the word 'noise' in the bottom right corner of this
% figure.  Compare with figure(4);   Note the weaker edges in
% the image are noisier with sigma = 1.

hFig = figure(1); clf;
image(noisyIm);
colormap(gray(256));
resizeImageFig(hFig, size(im), sclPix);

fprintf(2, 'Figures 2, 3  sigma = 2, 1 no noise.\n');
fprintf(2, 'Figures 4, 5  sigma = 2, 1 + sigma = %f noise.\n', nStdDev);
fprintf(2,...
  'The figure windows may all be on top of each other. Separate them.\n');
fprintf(2,'Notice the higher contrast edges are stable, %s\n', ...
           'despite the added noise.');

%% Conclude: 
%%   High contrast edgels are stable, despite the added noise.
%%   Lower constrast edgels are less stable.
fprintf(2,'Press any key to continue...\n'); pause;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Parkbench Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compare the runtime for the non-max suppression written using loops
%% with the runtime for an algorithm written without loops.
im = pgmRead('parkbench.pgm');
%% Display it
figure(1);clf;
showIm(im);
fprintf(2,'Press any key to continue...\n'); pause;

% Compute Canny edgels
sigma = 2.0;
hideSteps = TRUE;
slow = TRUE;  %% Uses loops in the non-max suppression
% Watch how long non-max supression takes...
[edgelIm0 nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, 4.0, ...
                                            hideSteps, slow);

%% Uses matrices instead of loops in the non-max suppression
% Watch how long non-max supression takes...
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, 4.0,...
                                            hideSteps, ~slow);

% Are the results the same?
% The following is zero if and only if the edgel images are identical
fprintf(2, 'Sanity check(xor edgels == 0?): %d\n', ...
   max(max(xor(edgelIm0, edgelIm))));

%% Fancy output: Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), sclPix);

%% Conclude: 
%%   Loops in uncompiled Matlab code are slow!
fprintf(2,'Press any key to continue...\n'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Zone Example %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A zone image is a common test case for image processing applications.
im = pgmRead('zone256.pgm');
stepHeight = 32;
noiseStdDev = 2.0;
disp(['zone amplitude: ' num2str(stepHeight) ...
      ', noiseStdDev: ' num2str(noiseStdDev)]);
im = 128 + (im/255.0-0.5) * stepHeight + randn(size(im)) * noiseStdDev;

%% Display it
figure(1);clf;
showIm(im);
title(sprintf('Noisy zone image (SNR = %f)', stepHeight/noiseStdDev));
fprintf(2,'Noisy zone image (SNR = %f/%f)\n', stepHeight, noiseStdDev);
fprintf(2,'Press any key to continue...\n'); pause;

%% Run Canny edges on this image, using a small sigma, and show the steps.
sigma = 1;
hideSteps = FALSE;  
[edgelIm nrmGradIm dirIm] = cannyEdgelDemo(im, sigma, 0.0, hideSteps);

%% Fancy output: Edgel image with a periodic colour map
enoughGradAmp = ~isnan(dirIm);
dirIm(isnan(dirIm)) = -0.1;
sclPix = round(max(10, 5000/max(size(im))))/10;
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* edgelIm + (~edgelIm)* 0.0);
resizeImageFig(hFig, size(dirIm), 2.0);

%% Note the edgels are noisy at the center of the disk.
%% This is due to small gradient amplitude in this region.
figure(1); clf;
showIm(nrmGradIm);

%% Also note the edgels are fragmented near the boundaries of the disk.
%% This is due to the low gradient amplitude, and also due to the
%% non-max suppresion.  Notice that the gradient direction image 
%% (without non-max suppression) is more coherent near the edge of the disk.
hFig = figure(1); clf;
colormap([0 0 0; hsv(64)]);
image( (ceil(dirIm*64/2.0)+1) .* enoughGradAmp + (~enoughGradAmp)* 0.0);
resizeImageFig(hFig, size(dirIm), 2.0);

%% Conclude: 
%%   The non-max suppression may be deleteing too many closely spaced
%%   (i.e. high frequency) edgels.
fprintf(2,'Press any key to end\n'); pause;



