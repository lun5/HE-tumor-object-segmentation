%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File: houghDemo.m  
%%% Matlab script file. Hough Transform for Lines Demo
%%%
%%% Compute Canny edgels using Gaussian derivative filter kernels,
%%% and non-max suppression.  Compute Hough transform of these edgels.
%%% It is left to the reader to play with bin sizes and thresholds to
%%% do image line detection.
%%%
%%% Dependencies
%%%   Directory: matlabVisTools
%%%      iseToolbox/pyrTools/
%%%          showIm.m  histo.m
%%%      utvisToolbox/file/
%%%          pgmRead.m resizeImageFig.m
%%%      utvisToolbox/edge/
%%%           cannyEdgels.m
%%%      localUtil/
%%%          doHough.m

%%% Things still to do:
%%%
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

lineRoot = [matlabVisRoot '/utvisToolbox/tutorials/lineTut']; 
addpath([lineRoot '/localUtil']);

% Check each required directories is on path (and some of the files)
which pgmRead

which cannyEdgels

%%%%%%%%%%%  Initialize random number generator %%%%%%%%%%%%%%%%%%%%%%%
% Random number generator seed:
seed = round(sum(1000*clock));
seed0 = seed;
% Restart the random number generators rand and randn
rand('state', seed0);
% We also need to start randn using the autmatically generated seedn:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

FALSE = (0==1);
TRUE = ~FALSE;
SUPERIMPOSE = TRUE;
SAVE = FALSE;
CROP_DEMO = FALSE; %  TRUE for cropped image from im9h.pgm

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
if FALSE  % Select one of the following images.
  
  im = pgmRead('three.pgm');
  im = pgmRead('microserf.pgm');
  im = pgmRead('wall.pgm');
  im = pgmRead('parkbench.pgm');
  
  %im = pgmRead('star_porch.pgm');
  
  fname = 'c:/images/boxes/boxes01-add';
  im = pgmRead([fname '.pgm']);

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

end

% Add some noise
sigmaGrey = 3.0;
im = im + sigmaGrey * randn(size(im));


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
cropBox = [ 1 1 size(im,2) size(im,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Pick an image region %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pts = [];
%%% For full image skip  to 'Set cropBox' below
if FALSE
  fprintf(2, ...
          'Mouse in a region from the edge map for further processing.\n');
  figure(hFig);
  pts = ginput(2); pts = round(pts)
elseif CROP_DEMO 
  pts = [39 88; 81 136];  % one example
end

%%%% Set cropBox
if (size(pts, 1) >=2)
  xBox = pts(:,1)';
  yBox = pts(:,2)';
  cropBox = [xBox(1) yBox(1) xBox(2) yBox(2)];
  if ((yBox(2) <= yBox(1)) || (xBox(2) <= xBox(1)))
    fprintf(2, 'Empty crop box: ');
    display(cropBox);
  else
    figure(hFig); hold on;
    plot([xBox(1) xBox(1) xBox(2) xBox(2) xBox(1)], ...
         [yBox(1) yBox(2) yBox(2) yBox(1) yBox(1)], ...
          'w');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Get edgels within selected region %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pEdgel, thetaEdgel] = getEdgelData(edgelIm, dirIm, cropBox);
% or use whole image: [p, theta] = getEdgelData(edgelIm, dirIm);

%%%%%%%
% We are going to ignore the edgel orientation theta in estimating
% the fitting line for now, as it is often noisier than the position
% information.  But this is definitely something we will consider
% including later.  
%%%%%%%

% Set the colormap for different orientations to be periodic.
hLine=figure(2); clf;
nOrient = 8;
colormap(hsv(nOrient));
if (size(cropBox,1) < 1)
  corners = [0 0; size(im)];
else
  corners = reshape(cropBox, 2,2)';
end
x0 = sum(corners, 1)'/2.0;
radius = sqrt(max(sum((corners - repmat(x0',2,1)).^2,1)));
% Set edgel labels according to their quantized orientation
lblEdgel = round(thetaEdgel/pi*nOrient);
% Wrap the labels according to the periodic range of orientations [0, pi].
lblEdgel = mod(lblEdgel, nOrient) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Plot edgels within selected region %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To avoid plotPoints being too slow, limit the number
% of edgels to be plotted to be at most MAX_PLOT
MAX_PLOT = 5000;
if size(pEdgel,1) > MAX_PLOT
  % Randomly sample at MAX_PLOT edgel positions
  idq = randperm(size(pEdgel,1));
  qEdgel = pEdgel(idq(1:MAX_PLOT), :);
  qLblEdgel = lblEdgel(idq(1:MAX_PLOT));
else
  qEdgel = pEdgel;
  qLblEdgel = lblEdgel;
end

figure(2); clf;
plotPoints(qEdgel, figure(2), qLblEdgel); 


%
% The code below sets the number of bins to be used for theta and r,
% then calls the Hough transform function
%

ntBin = 60;
rRange = [-radius, radius];
dr = 2;
nrBin = max((rRange(2)-rRange(1))/dr, 1); 
[H tBins rBins] = doHough(pEdgel, x0, ntBin, nrBin, rRange, 0);
 
%%%%%  Explore with point and click

figure(1); close 1; figure(1);
showIm(H); resizeImageFig(figure(1), size(H), 3);

% You may wish to blow up a portion of figure(1) now, near some
% of the local maxima.
% Also, if the Hough transform image looks too dark, try plotting the
% log of the Hough transform, i.e. showIm(log(H+k));
% k is a small constant (.1 - 10) used to avoid computing log of 0

% Take a look at a plot of the Hough Transform surface, if the transform
% isn't too large...
if prod(size(H)) < 200^2
  figure(5); clf; surf(H); colormap(jet); shading interp;
end
% Here you should be able to see the peaks and valleys of the
% function.  If we find the big peaks in this plot, do we get all
% the lines in the original image, and only the salient lines?  We
% investigate this below.
fprintf(2,'Press any key to continue...\n'); pause;


% Use the code below to see the lines in image space that correspond to particular
% bins in the Hough transform image, try finding the points in the Hough transform
% that correspond to image lines. Also, notice what happens if you click on regions
% that are not local maxima. Are there any points that look like local maxima
% yet do not correspond to a line in the image?

figure(2); clf; plotPoints(qEdgel,figure(2), qLblEdgel);
%%% Point and click ten times... or use CTRL_C to break out
for kLine=1:10
  figure(1);
  fprintf(2, ' Mouse in a Hough peak in Fig. 1 ...');
  tr = ginput(1)+0.5; %% Shift pixel borders to integer positions
  if isempty(tr)
    fprintf(2, ' No moused point provided.  Break.\n');
    break;
  end
  tr = floor(tr);
  tr(1) = max(tr(1),1); tr(1) = min(tr(1), size(H,2));
  tr(2) = max(tr(2),1); tr(2) = min(tr(2), size(H,1));

  tWin = tr(1)+ (-2:2);
  tWin = tWin(tWin >0 & tWin < size(H,2));
  rWin = tr(2)+ (-2:2);
  rWin = rWin(rWin >0 & rWin < size(H,1));
  H(rWin, tWin)

  lineParams = [tBins(tr(1)) rBins(tr(2))];
  nrPar = [cos(lineParams(:,1)), sin(lineParams(:,1)), ...
           zeros(size(lineParams,1),1)];
  nrPar(:,3) = lineParams(:,2) - nrPar(:,1:2) * x0;

  h =figure(2);  hold on;
  drawLines(nrPar, cropBox);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now, see what happens when we perform non-maximal suppresion on the Hough
% transform array

nonMaxRadius = 3;
[Hloc tBins rBins] = doHough(pEdgel, x0, ntBin, nrBin, rRange, ...
                          nonMaxRadius);
figure(4); clf; showIm(Hloc);

%%%%%%%%%%%%%% Hough Histograms
% Recall the Gradient Magnitude histogram we used to determine the threshold
% for canny edgels. 
%%%%%%%%%%%%%%

meanVotes = size(pEdgel,1)*length(tBins)/(length(tBins)*length(rBins))

[n x] = histo(H(H>0), 64);
figure(3); clf; plot(x,n, 'b');

[n x] = histo(Hloc(Hloc>0), 64);
figure(3); hold on; plot(x,n, 'r');


%%%%%%%%%%%%%% Hough Threshold
% Find the point in the Hough Transform that is close to the top 95
% percentile
%%%%%%%%%%%%%%

Hprc = prctile(H(H>0), 95)
Hmin = Hprc;   %%%%%%%%%%%% TRY CHANGING THIS threshold Hmin
q = find(Hloc>Hmin);

%%%%%%%%%%%%%
% Finally, plot the lines that correspond to all bins in the Hough Transform
% that are above the specified threshold
%%%%%%%%%%%%%

lineParams = [];
if length(q)>0
  fprintf(2, ' %d lines above threshold %f\n', length(q), Hmin);
  qt = floor((q-1)/size(H,1))+1;
  qr = q - (qt-1)*size(H,1);
  
  lineParams = [tBins(qt)  (rBins(qr)+dr/2)];

  nrPar = [cos(lineParams(:,1)), sin(lineParams(:,1)), ...
           zeros(size(lineParams,1),1)];
  nrPar(:,3) = lineParams(:,2) - nrPar(:,1:2) * x0;

  h =figure(4); clf; plotPoints(qEdgel, h); hold on;
  drawLines(nrPar, cropBox);

else
  fprintf(2, ' No lines above threshold %f\n', Hmin);
end  

%%% Variations: 
% Try setting CROP_DEMO = false 

