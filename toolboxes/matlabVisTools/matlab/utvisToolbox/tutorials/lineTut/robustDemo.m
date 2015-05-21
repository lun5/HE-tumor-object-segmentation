%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File: robustDemo.m  
%%% Matlab script file. Robust M-estimation for Lines Demo
%%%
%%% Compute Canny edgels using Gaussian derivative filter kernels,
%%% and non-max suppression. 
%%% For demonstration purposes (ONLY) we compute the robust objective 
%%% function for estimating (infinite) lines, and display it.
%%% A point and click interface is given for the iterative least squares
%%% algorithm discussed in class.
%%%
%%% Dependencies
%%%   Directory: matlabVisTools
%%%      iseToolbox/pyrTools/
%%%          showIm.m  histo.m
%%%      utvisToolbox/file/
%%%          pgmRead.m resizeImageFig.m
%%%      utvisToolbox/edge/
%%%           cannyEdgels.m computeOrientTensor.m
%%%      ./localUtil/

%%% ADJ, Oct. '03.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;

global matlabVisRoot

% We need to ensure the path is set for the iseToolbox.
if isempty(matlabVisRoot)
  dir = pwd;
  cd ~jepson/pub/matlab   % CHANGE THIS
  startup;
  cd(dir);
end

lineRoot = [matlabVisRoot '/utvisToolbox/tutorials/lineTut']; 
addpath([lineRoot '/localUtil']);

% Check each required directories is on path (and some of the files)
which pgmRead
which cannyEdgels

%%%%%%%%%%  Initialize random number generator %%%%%%%%%%%%%%%%%%%%%%%
% Random number generator seed:
seed = round(sum(1000*clock));
seed0 = seed;
% Restart the random number generators rand and randn
rand('state', seed0);
% We also need to start randn using the autmatically generated seedn:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

global FALSE;
global TRUE; 
FALSE = (0==1);
TRUE = ~FALSE;
SUPERIMPOSE = TRUE;
SAVE = FALSE;
CROP_DEMO = TRUE; %TRUE for initial demo using cropped image from im9h.pgm

% For display, scale small images up to be about 500 pixels
% in the largest dimension.  Larger images won't be scaled down.
figSz = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Blocks-World Example %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = pgmRead('im9h.pgm');
sclPix = round(max(10, (figSz * 10)/max(size(im))))/10;
% Display it
figure(1);close 1; figure(1); clf;
showIm(im);
fprintf(2,'Press any key to continue...\n'); pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Other image files  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FALSE  % Select one of the following images.
  
  im = pgmRead('three.pgm');
  im = pgmRead('microserf.pgm');
  im = pgmRead('wall.pgm');
  im = pgmRead('parkbench.pgm');
  
  % ./localUtil/
  im = pgmRead('star_porch.pgm');
  im = pgmRead('boxes01-add.pgm');

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
sigmaGrey = 0.0;
im = im + sigmaGrey * randn(size(im));


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

% Crop whole image by default
cropBox = [ 0 0 size(im,2) size(im,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Pick an image region %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pts = [];
% For full image skip  to 'Set cropBox' below
if FALSE
  fprintf(2, ...
          'Mouse in a region from the edge map for further processing.\n');
  figure(hFig);
  pts = ginput; pts = round(pts)
elseif CROP_DEMO 
  pts = [39 88; 81 136];  % one example
end

%%% Set cropBox
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
%%%%%%%%%%%%%% Get edgels within selected region %%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%% Plot edgels within selected region %%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Robust Objective Function %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We wish to generate an image of the objective function, darker
% values will have a lower (more optimal) value.

sigmaRho = 1;
ntBin = 60;
rRange = [-radius, radius];
dr = 1.5 * sigmaRho;
nrBin = max(ceil((rRange(2)-rRange(1))/dr), 1); 
[E tBins rBins] = sampleRobustObj(pEdgel, x0, sigmaRho,...
                                  ntBin, nrBin, rRange, 0);

% Display objective function
figure(1); close 1; figure(1);
showIm(E); resizeImageFig(figure(1), size(E), 3);
% Here the horizontal axis corresponds to the angle
% of the normal, theta, in radians, and it goes from 0 to pi.
% The vertical axis is the perpendicular distance of the line from the
% center point x0.  The range of the veritcal axis is given by
% [-radius, +radius] (increasing downwards) where radias is the
% radius of the image in Figure 2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Explore Objective Function %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redraw image of objective function 
figure(1); close 1; figure(1);
showIm(E); resizeImageFig(figure(1), size(E), 3);
% You may wish to blow up a portion of figure(1) now, near some
% of the local minima.
fprintf(2,'Press any key to continue...\n'); pause;

% Take a look at a plot of the Hough Transform surface, if the transform
% isn't too large...
if prod(size(E)) < 200^2
  figure(5); clf; surf(E); colormap(jet); shading interp;
  fprintf(2,'Press any key to continue...\n'); pause;
end
% Here you can easily see the plateaus and valleys of the error
% function. See what happens when you plot the error surface
% for the complete image.

% Use the left mouse button to select points in the image
% of the objective function.  These points correspond to particular
% values of theta and r.  The corresponding line:
%   n(theta) dot (x - x0) + r = 0
% Is drawn in figure(2), the edgel plot.
% If figure(2) has been corrupted, redraw it:
%figure(2); clf; plotPoints(qEdgel,figure(2), qLblEdgel);
% Use CTRL_C to break out of loop early if you wish.
for kLine=1:10
  
  % Wait for a mouse click.
  figure(1);  fprintf(2, ' Mouse in a point in Fig. 1 ...');
  tr = ginput(1)+0.5; %% Shift pixel borders to integer positions
  if isempty(tr)
    fprintf(2, ' No moused point provided.  Break.\n');
    break;
  end
  tr = floor(tr);
  tr(1) = max(tr(1),1); tr(1) = min(tr(1), size(E,2));
  tr(2) = max(tr(2),1); tr(2) = min(tr(2), size(E,1));

  % Dump out objective function in 5x5 neighbourhood 
  % centered on clicked point
  tWin = tr(1)+ (-2:2);
  tWin = tWin(tWin >0 & tWin < size(E,2));
  rWin = tr(2)+ (-2:2);
  rWin = rWin(rWin >0 & rWin < size(E,1));
  E(rWin, tWin)

  % Get the line parameters for clicked point, and
  % draw the line in Figure 2.
  lineParams = [tBins(tr(1)) rBins(tr(2))];
  nrPar = [cos(lineParams(:,1)), sin(lineParams(:,1)), ...
           zeros(size(lineParams,1),1)];
  nrPar(:,3) = lineParams(:,2) - nrPar(:,1:2) * x0;

  h =figure(2);  hold on;
  drawLines(nrPar, cropBox);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Local Minima of Objective Function %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nonMaxRadius = 3;
% Find local minima within a radius of nonMaxRadius...
[Eloc tBins rBins] = sampleRobustObj(pEdgel, x0, sigmaRho,...
                            ntBin, nrBin, rRange, nonMaxRadius);
% Plot the local min image
figure(4); clf; 
showIm(Eloc); resizeImageFig(figure(4), size(E), 3);

%%%%%%%%%%%%% Histogram of Objective Function %%%%%%%%%%%%%%%
% The maximum of the objective function is equal
% to the maximum of the robust estimator rho, which is one,
% times the number of edgels size(pEdgel,1). 
[n x] = histo(E(:), 64);
figure(3); clf; plot(x,n, 'b');
% There could be a large peak in the neighbourhood of the
% maximum possible value, which causes the plot to be poorly
% scaled.  Avoid this peak as follows...
[n x] = histo(E(E<size(pEdgel,1)-2), 64);
figure(3); clf; plot(x,n, 'b');

[n x] = histo(Eloc(Eloc<size(pEdgel,1)-2), 64);
figure(3); hold on; plot(x,n, 'r');
title('Histogram of objective fcn E (b), local min (r)');

fprintf(2,'Press any key to continue...\n'); pause;
figure(3); close 3;


%%%%%%%%%%%%% Objective Function Threshold
Eprc = prctile(E(E<size(pEdgel,1)-2), 60)
Emax = Eprc;  %% On more complex images, play with this threshold

q = find(Eloc<Emax);
fprintf(2, 'Selected %d local min.\n', length(q));

%%%%%%%%%%%%% Draw corresponding lines, if any
lineParams = [];
if length(q)>0
  fprintf(2, ' %d lines below threshold %f\n', length(q), Emax);
  qt = floor((q-1)/size(E,1))+1;
  qr = q - (qt-1)*size(E,1);
  
  lineParams = [tBins(qt)  (rBins(qr)+dr/2)];

  nrPar = [cos(lineParams(:,1)), sin(lineParams(:,1)), ...
           zeros(size(lineParams,1),1)];
  nrPar(:,3) = lineParams(:,2) - nrPar(:,1:2) * x0;

  figure(4); close 4; h=figure(4); clf; plotPoints(qEdgel, h); hold on;
  drawLines(nrPar, cropBox);

else
  fprintf(2, ' No lines below threshold %f\n', Emax);
end 
fprintf(2,'Press any key to continue...\n'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Local Minima at a Higher Resolution %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It is instructive to look at the objective function at
% a higher resolution.  Since it is a continuous function, we
% only need to increase the number of bins in theta and r
% and replot.  Because of the time it takes to fully sample
% a high resolution image, we will only do this on the
% crop-demo example (run the above code with CROP_DEMO=TRUE)
if CROP_DEMO
  % Increase the resolution by 8.
  ntBin = 60 * 4;
  dr = 1.5 * sigmaRho/ 8;
  nrBin = max(ceil((rRange(2)-rRange(1))/dr), 1); 
  [E tBins rBins] = sampleRobustObj(pEdgel, x0, sigmaRho,...
                                    ntBin, nrBin, rRange, 0);
  sclFig = round(figSz/max(size(E))); sclFig = max(sclFig, 1);
  figure(1); close 1; figure(1);
  showIm(E); resizeImageFig(figure(1), size(E), sclFig);
  % Beautiful!....
  
  % Find local minima within a radius of nonMaxRadius...
  nonMaxRadius = 1;
  [Eloc tBins rBins] = sampleRobustObj(pEdgel, x0, sigmaRho,...
                                       ntBin, nrBin, rRange, nonMaxRadius);
  % Plot the local min image
  figure(4); clf; 
  showIm(Eloc); resizeImageFig(figure(4), size(E), sclFig);
  % Look carefully... there are many local minima.  For
  % Better contrast...
  Etmp = Eloc;
  Etmp(Etmp == size(pEdgel,1)) = -1;
  figure(4); clf; 
  showIm(Etmp); resizeImageFig(figure(4), size(E), sclFig);
  
  % Use a high threshold to show the lines for all local minima.
  Emax = size(pEdgel,1)-1;  

  q = find(Eloc<Emax);
  fprintf(2, 'Selected %d local min.\n', length(q));

  %%%%%%%%%%%%% Draw corresponding lines, if any
  lineParams = [];
  if length(q)>0
    fprintf(2, ' %d lines below threshold %f\n', length(q), Emax);
    qt = floor((q-1)/size(E,1))+1;
    qr = q - (qt-1)*size(E,1);
    
    lineParams = [tBins(qt)  (rBins(qr)+dr/2)];

    nrPar = [cos(lineParams(:,1)), sin(lineParams(:,1)), ...
             zeros(size(lineParams,1),1)];
    nrPar(:,3) = lineParams(:,2) - nrPar(:,1:2) * x0;

    figure(4); close 4; h=figure(4); clf; plotPoints(qEdgel, h); hold on;
    drawLines(nrPar, cropBox);

  else
    fprintf(2, ' No lines below threshold %f\n', Emax);
  end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Important Observation:
  % The point here is that the objective function has a few
  % deep valleys, corresponding to image lines, along with
  % several plateaus.  These plateaus have weak local minima,
  % presumably due to small local alignments of edgels (often
  % involving distant edgels), which cause ripples in the values
  % on the plateaus... and hence, occasionally, weak local
  % minima.  The importance of this observation is that if
  % we use a local hill-descent algorithm, we can end up
  % getting stuck in these local minima.  (Note that these
  % are probably not ALL the local minima, either, since
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % Clean up.
  clear Etmp;
  % Reset original objective function parameters
  ntBin = 60;
  dr = 1.5 * sigmaRho;
  nrBin = max(ceil((rRange(2)-rRange(1))/dr), 1); 
  nonMaxRadius = 3;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Robust M-estimation Algorithm %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Explore with point and click
[E tBins rBins] = sampleRobustObj(pEdgel, x0, sigmaRho,...
                                  ntBin, nrBin, rRange, 0);


nIts = 200;  % Maximum number of iterations
figure(4); close 4;

% Do initial plots for point and click
figure(1); close 1; figure(1);
sclFig = round(figSz/max(size(E))); sclFig = max(sclFig, 1);
showIm(E); resizeImageFig(figure(1), size(E), sclFig);
h=figure(4); clf; plotPoints(qEdgel, h); hold on;

% Check out the convergence of the robust iteration algorithm.
for kLine=1:10
  figure(1);
  fprintf(2, 'Select initial point in Fig. 1 ...\n');
  tr = ginput(1)+0.5; %% Shift pixel borders to integer positions
  if isempty(tr)
    fprintf(2, ' No moused point provided.  Break.\n');
    break;
  end
  % Compute nearest pixel center
  tr = floor(tr);
  tr(1) = max(tr(1),1); tr(1) = min(tr(1), size(E,2));
  tr(2) = max(tr(2),1); tr(2) = min(tr(2), size(E,1));
  % Plot start position as a green circle.
  figure(1); hold on; plot(tr(1), tr(2), 'og');

  % Compute corresponding line parameters and call robustIteration alg.
  lineParams = [tBins(tr(1)) rBins(tr(2))];
  nrPar = [cos(lineParams(:,1)), sin(lineParams(:,1)), ...
           zeros(size(lineParams,1),1)];
  nrPar(:,3) = lineParams(:,2) - nrPar(:,1:2) * x0;

  n0 = nrPar(1,1:2)'; r0 =  nrPar(1,3);
  
  [n, r, converged, nVals, rVals, sumWs ] = ...
      robustIteration(n0, r0, pEdgel, sigmaRho, nIts);
  
  % Report results from robustIteration alg.
  if converged
    fprintf(2,' Converged in %d iterations\n', size(nVals,2));
  else
    fprintf(2,' Has not yet converged after %d iterations\n', ...
            size(nVals,2));
  end
  figure(2); clf; plot(sumWs);
  xlabel('Iteration'); ylabel('Sum of weights');
  
  % Display iterations of robustIteration alg. on objective
  % function plot.
  tE = atan2(nVals(2,:), nVals(1,:));
  rE = rVals;
  idx = tE >= pi;
  if any(idx)
    tE(idx) = tE(idx) - pi;
    rE(idx) = -rE(idx);
  end
  idx = tE < 0;
  if any(idx)
    tE(idx) = tE(idx) + pi;
    rE(idx) = -rE(idx);
  end
  nVals = [cos(tE); sin(tE)];
  rE = rE + x0' * nVals; 
  dt = tBins(2) - tBins(1);
  dr = rBins(2) - rBins(1);
  tr = [((tE - tBins(1))/dt +1) ; ((rE-rBins(1))/dr + 1)];
  figure(1);
  if size(tr,2)>1
    plot(tr(1,2:end), tr(2,2:end), 'ob');
  end
  if converged
    plot(tr(1,end), tr(2,end), '*r');
  end
  
  % Display iterations of robustIteration alg. on plot of
  % edgel points.
  nrPar = [cos(tE(1,:))', sin(tE(1,:))', ...
           zeros(size(tE,2),1)];
  nrPar(:,3) = rE(:) - nrPar(:, 1:2) * x0;

  h=figure(4); clf; plotPoints(qEdgel, h); hold on;
  drawLines(nrPar, cropBox);
  fprintf(2,'Press any key to continue...\n'); pause;

end
fprintf(2, 'Done convergence trials.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations:
% The algorithm appears to converge given an intial guess
% close to a strong local minima.  However, given more distant
% initial guesses, it can get hung up on weak local minima which
% do not correspond to any salient line in the edgel data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rerun this with more cluttered images.  For example, 
% turn off the cropping (CROP_DEMO=FALSE) and run on the
% whole im9h.pgm image.  Or use some of the other images
% (cropped or not) provided.
 
