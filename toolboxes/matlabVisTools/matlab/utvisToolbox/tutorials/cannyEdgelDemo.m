function [edgelIm, nrmGradIm, dirIm] = cannyEdgelDemo(im, sigma, ...
                                          minStrength, quiet, slow)
% [edgelIm nrmIm dirIm] = cannyEdgels(IM, SIGMA, MINSTRENGTH, QUIET, SLOW)
%
% Compute the Canny edgels using the 1st derivatives of a Gaussian filter,
% and display the intermediate results.  This does not perform 
% hysteresis thresholding.
%
% PARAMS:  
%  IM the 2D gray level input image
%  SIGMA (optional, default = 1.0) the std dev of the Gaussian filter.
%  MINSTRENGTH (optional, default = moused in value)  The threshold
%    value for the minimum gradient magnitude.  
%    If minStrength <= 0.0, then a histogram of the gradient magnitudes
%    is plotted, and the user is asked to mouse in a value.
%  QUIET (optional, default = 1 (quiet)) A zero value causes several
%        intermediate results to be plotted.
%  SLOW (optional, default = 0 (not slow)) a nonzero value causes
%       the non-max suppression to be done using a slow algorithm
%       which loops through the image pixels.  Otherwise a faster
%       matrix based algorithm is used. 
%
% OUTPUT:  
%  edgelIm  binary image marking edgel locations.
%  nrmGradIm  image of norm of gradient image (NOT thresholded)
%  dirIm    direction of gradient image (in terms of theta/pi), so
%           directions are periodic in the interval (0, 2]. 
%           dirIm is NaN at pixels for which nrmGradIm < minStrength.  

% ADJ 9/01.

%%%%%%%%%%%%%%% Check parameters  %%%%%%%%%%%%%%%%%%%%%%%

%%% Fill in default parameter values, and correct bogus parameter values.
if ~exist('sigma', 'var')  %% sigma is the std dev for the Gaussian filters
  sigma = 1.0;
elseif (sigma < 0.5/3)
  fprintf(2, 'User specified sigma = %e is too small.\n', sigma); 
  sigma = 0.5000001/3;
  fprintf(2, 'Using sigma = %f instead\n', sigma);
end

if ~exist('minStrength', 'var')  
  useMouse = 1;            %% indicates that the minStrength threshold
                           %% is to be moused in.
elseif minStrength<=0.0
  useMouse = 1;
else
  useMouse = 0;            %% minStrength provided as an argument.
  if (minStrength < eps)
    fprintf(2,'Invalid edge strength tolerance %e\n', minStrength); 
    minStrength = eps;
    fprintf(2,'Using edge strength tolerance = %e\n', minStrength);
  end
end

if ~exist('quiet', 'var')
  quiet = 1;              %% do not display intermediate results.
end

if ~exist('slow', 'var')
  slow = 0;               %% use the quicker non-max suppression code.
end 


% For display, scale small images up to be about 500 pixels
% in the largest dimension.  Larger images won't be scaled down.
figSz = 500;
sclPix = round(max(10, (figSz * 10)/max(size(im))))/10;

%%%%%%%%%%%%%%% Build Filter Kernels %%%%%%%%%%%%%%%%%%%%%%%

%%% Build Gaussian filter masks, along with derivatives.
%%% The second derivative is not used yet.
sigmaSqr = sigma*sigma;
gFiltSize = 2 * round(3.0 * sigma) + 1;
x = [1:gFiltSize] - round((gFiltSize+1)/2);
gFilt = exp(- x .* x / (2.0*sigmaSqr));
gFilt = gFilt/ sum(gFilt(:));
gxFilt = - (x / sigmaSqr) .* gFilt;
gxxFilt = ((x / sigmaSqr).^2 - 1.0/sigmaSqr) .* gFilt;

if ~quiet
  %%% Plot the 1D filters
  figure(1); close; figure(1);
  plot(x, gFilt, 'b');
  hold on;
  plot(x, gxFilt, 'g');
  plot(x, gxxFilt, 'r');
  title('Gaussian Filter (b), First and Second Derivatives (g, r)');
  xlabel('x');
  ylabel('Filter Tap: filt(x)');
  hold off;
  fprintf(2, 'Displaying Gaussian derivative filters for sigma = %f\n', ...
          sigma);
  fprintf(2,'Press any key to continue...\n'); pause;
end

%%%%%%%%%%%%%%% Do separable convolutions %%%%%%%%%%%%%%%%%%%

gradIm = zeros([size(im), 2]);
gradIm(:,:,1) = rconv2sep(im, gxFilt, gFilt);
gradIm(:,:,2) = rconv2sep(im, gFilt, gxFilt);

%%Type 'help rconv2sep'  for more information.

if ~quiet
  figure(1); clf
  subplot(1,2,1);
  showIm(gradIm(:,:,1));
  title('G_X * I');
  subplot(1,2,2);
  showIm(gradIm(:,:,2));
  title('G_Y * I');
  fprintf(2,'Displaying x,y gradient images\n');
  fprintf(2,'Press any key to continue...\n'); pause;
end

%%%%%%%%%%%%%%% Gradient Magnitude Image %%%%%%%%%%%%%%%%%%%
 
nrmGradIm = sum(gradIm .* gradIm, 3).^0.5;

if ~quiet
  figure(1); clf;
  showIm(nrmGradIm);
  title('Image Gradient Magnitude');
  fprintf(2, 'Displaying gradient magnitude image\n');
  fprintf(2,'Press any key to continue...\n'); pause;
end

%%%%%%%%%%%%%%% Select Magnitude Threshold %%%%%%%%%%%%%%%%%%%
if useMouse
   %% Build a histogram of the gradient magnitudes.
  [hist bin] = histo(nrmGradIm, 64);

  %% Plot the histogram
  figure(1); clf;
  plot(bin, log10(hist + 1));
  title('Log histogram of amplitude');
  xlabel('Gradient Amplitude')
  ylabel('Log Frequency')
  axis tight;
  drawnow;
 
  %% Select a threshold value on this histogram.
  %% In many cases the histogram has a peak of responses for small
  %% amplitudes, and then a linearly flat or decaying tail.  If your
  %% image has this stucture, select a point just beyond this peak, 
  %% at the beginning of the tail.
  fprintf(2, 'Mouse in the low-amplitude threshold (See comments\n');
  fprintf(2, 'in cannyEdgelDemo.m for details)...\n');
  thres = ginput;
  %% Use the last value moused in (if any)
  if (size(thres,1) < 1)
    minStrength = 0.0;  % Nothing moused in.
  else
    minStrength = thres(size(thres,1), 1);
  end
  if (minStrength < eps)
    minStrength = eps;
  end
  %% Display selected threshold.
  hold on;
  ax = axis;
  plot([minStrength minStrength], ax(3:4), 'r');
  hold off;
  fprintf(2, ' Minimum edge strength tolerance: %f\n', minStrength);
  fprintf(2,'Press any key to continue...\n'); pause;
end;

%%%%%%%%%%%%%%% Enforce Magnitude Threshold %%%%%%%%%%%%%%%%%%%

dirIm = zeros(size(im));
enoughGradAmp = (nrmGradIm >= minStrength);

if ~quiet
  figure(1); clf
  subplot(1,2,1);
  showIm(nrmGradIm);
  title('Norm of Image Gradient');
  subplot(1,2,2);
  showIm(double(enoughGradAmp));
  title('Thresholded Region');
  fprintf(2,'Displaying norm of gradient, and thresholded image region.\n');
  fprintf(2,'Press any key to continue...\n'); pause;
end

%%%%%%%%%%%%%%% Compute Direction Image %%%%%%%%%%%%%%%%%%%
dxIm = gradIm(:,:,1);
dyIm = gradIm(:,:,2);
dirIm(enoughGradAmp) = atan2(dyIm(enoughGradAmp), dxIm(enoughGradAmp))/pi;
%% Note this value of dirIm is periodic in the range [-1, 1]

if ~quiet
  %%% Display direction image
  figure(1); clf;
  blankVal = -1.2;
  showIm(dirIm .* enoughGradAmp + (~enoughGradAmp)*blankVal);
  title('Gradient Direction: [-\pi,\pi]');
  fprintf(2, 'Displaying gradient direction image\n');
  drawnow; 
  fprintf(2,' Branch cut at 180 degrees\n');
  fprintf(2,'Press any key to continue...\n'); pause;

  %% You may notice a discontinuity in the dirIm for nearly
  %% vertical edges.  The value of dirIm hops between -1 and 1.
  %% This discontinuity is an artifact of plotting the
  %% direction as a gray level.  The discontinuity
  %% is at the so-called "branch cut" for the angle.

  %% For example we can move the branch cut to 45 degrees 
  %% (where the corresponding edges look like "dark / bright")
  cut = 45/180;
  while (sum(sum(dirIm < cut)) > 0)
    dirIm(dirIm<cut) = dirIm(dirIm<cut) + 2.0;
  end
  while (sum(sum(dirIm >= cut+2.0)) > 0)
    dirIm(dirIm>=cut+2.0) = dirIm(dirIm>=cut+2.0) - 2.0;
  end
  figure(1); clf;
  showIm(dirIm .* enoughGradAmp + (~enoughGradAmp)*blankVal);
  title(sprintf('Gradient Direction: [%5.2f\\pi,%5.2f\\pi]',cut, cut+2));
  drawnow; 
  fprintf(2,' Branch cut at %f degrees\n', cut*180);
  fprintf(2,'Press any key to continue...\n'); pause;
  %% Now the vertical edges appear smooth, but not the ones near 45 degrees.

  %% Branch cut at 0 degrees (edges look like "dark | bright")
  cut = 0;
  while (sum(sum(dirIm > cut+2.0)) > 0)
    dirIm(dirIm>cut+2.0) = dirIm(dirIm>cut+2.0) - 2.0;
  end
  while (sum(sum(dirIm <= cut)) > 0)
    dirIm(dirIm<=cut) = dirIm(dirIm<=cut) + 2.0;
  end
  figure(1); clf;
  showIm(dirIm .* enoughGradAmp + (~enoughGradAmp)*blankVal);
  title(sprintf('Gradient Direction: [%5.2f\\pi,%5.2f\\pi]',cut, cut+2));
  drawnow; 
  fprintf(2,' Branch cut at %f degrees\n', cut*180);
  fprintf(2,'Press any key to continue...\n'); pause;

  %% Better yet, use a colour map that is periodic.
  hFig = figure(1); clf;
  colormap([0 0 0; hsv(64)]);
  image( (ceil(dirIm*64/2.0)+1) .* enoughGradAmp + (~enoughGradAmp)* 0.0);
  resizeImageFig(hFig, size(dirIm), sclPix);
  fprintf(2,' Direction image with a periodic colour map\n');
  fprintf(2,'Press any key to continue...\n'); pause;

end

%% Branch cut at 0 degrees (edges look like "dark | bright")
cut = 0;
while (sum(sum(dirIm > cut+2.0)) > 0)
  dirIm(dirIm>cut+2.0) = dirIm(dirIm>cut+2.0) - 2.0;
end
while (sum(sum(dirIm <= cut)) > 0)
  dirIm(dirIm<=cut) = dirIm(dirIm<=cut) + 2.0;
end
%%% dirIm now in the range (0, 2.0]

%% Insert NaN in regions for which the gradient is too small.
dirIm(~enoughGradAmp) = NaN;


%%%%%%%%%%%%%%% Do Non-max Edgel Suppression %%%%%%%%%%%%%%%%
edgelIm = zeros(size(im));
%% Set alpha to be the distance from the origin to the
%% point that is one pixel to the right, and pi/8 degrees off horizontal.
theta = pi/8;
alpha = 0.5/sin(theta);
%% We scale the edge normal by this amount, and then round it.
%% This causes the rounded value to be the closest direction
%% vector on the discrete grid for a 3x3 neighbourhood.
if ~quiet
  %%% Display direction image
  figure(1); clf;
  t = [0:0.01:1] * 2 * pi;
  hold on;
  for k = -1:1
    plot([-1.5 1.5], [k k], 'k');
    plot([k k], [-1.5 1.5], 'k');
  end
  for k = -0.5:0.5
    plot([-1.5 1.5], [k k], 'k:');
    plot([k k], [-1.5 1.5], 'k:');
  end
  axis([-1.5 1.5 -1.5 1.5]);
  ha = get(gcf, 'CurrentAxes');
  plot(alpha*cos(t), alpha*sin(t), 'r');
  samplePts = [-1 -1; -1 0; -1 1; ...
               0 -1;  0 0;  0 1; ...
               1 -1;  1 0;  1 1];
  scatter(samplePts(:,1), samplePts(:,2), 'ko');
  for k = 1:8
    plot([0 alpha*cos((2*k+1)*theta)], [0 alpha*sin((2*k+1)*theta)], 'g');
  end
  hold off;
  title('Rounding Gradient Directions to Grid');
  xlabel('Image x coord');
  ylabel('Image y coord');
  fprintf(2,'Press any key to continue...\n'); pause;
end

%%
%% Loop through the image, checking for larger amplitude
%% gradient values in the neighbourhood of each pixel, in
%% the direction of the gradient (rounded to the nearest pixel).
%% If a larger value is found in either direction, then the
%% current pixel is not an edgel.

fprintf(2,'Begin non-maximal edgel supression: ');
if slow
  fprintf(2,' Loop version...'); tic;
  %%%%%% Slow, but 'easy to read' non-max suppression implementation
  for x=1:size(im,2);
    for y=1:size(im,1);
      if enoughGradAmp(y,x)
        ok = 1;
        q = round(alpha * [dyIm(y,x) dxIm(y,x)] /nrmGradIm(y,x));
        p = [y x] + q;
        if ((min(p) > 0) & (p(1) <= size(im,1)) & (p(2) <= size(im,2)))
          if ((enoughGradAmp(p(1), p(2))) & ...
              (nrmGradIm(p(1), p(2)) > nrmGradIm(y,x)))
            ok = 0;
          end
        end
        if ok
          p = [y x] - q;
          if ((min(p) > 0) & (p(1) <= size(im,1)) & (p(2) <= size(im,2)))
            if ((enoughGradAmp(p(1),p(2))) & ...
                (nrmGradIm(p(1),p(2)) > nrmGradIm(y,x)))
              ok = 0;
            end
          end
        end
        edgelIm(y,x) = ok;
      end
    end
  end

else

  fprintf(2,' Matrix ops version...'); tic;
  %%%%%% Matrix version of non-max suppression implementation

  %% Make the norm of the gradient artificially large at
  %% pixels that are below the minStrength threshold.
  maxGrad = max(nrmGradIm(:));
  nrm = nrmGradIm .* enoughGradAmp + (~enoughGradAmp) * (10.0 * maxGrad);
  
  %% Compute the shifts in the x and y directions.
  qx = round(alpha * dxIm ./ nrm);
  qy = round(alpha * dyIm ./ nrm);

  %% Interpolate the norm of the gradient image at the
  %% pixels given by the shifts.
  [x y] = meshgrid(1:size(im,2), 1:size(im,1));
  ampNear = interp2(x, y, nrmGradIm, x+qx, y+qy, 'nearest');
  ampNear(isnan(ampNear)) = 0.0;
  
  %% Do non-max supression in this shift direction.
  ok = enoughGradAmp;
  ok = ok & (ampNear < nrmGradIm);

  %% Do non-max supression in the opposite shift direction.
  ampNear = interp2(x, y, nrmGradIm, x-qx, y-qy, 'nearest');
  ampNear(isnan(ampNear)) = 0.0;
  edgelIm = ok & (ampNear < nrmGradIm);
end
fprintf(2,'done\n'); toc;
 
if ~quiet  
  %%%%%%%%%%%%%%%%%% Display edgel image %%%%%%%%%%%%%%%%%%
  hFig = figure(1); clf;
  image(edgelIm*255);
  colormap(gray(256));
  resizeImageFig(hFig, size(edgelIm), sclPix);
  fprintf(2,' Displaying binary edgel image\n');
  fprintf(2,'Press any key to continue...\n'); pause;

  overlay = 0.8 * im;
  overlay(edgelIm>0) = 255;
  image(overlay);
  resizeImageFig(hFig, size(edgelIm), sclPix);
  fprintf(2,' Overlaying edgels on original image\n');
  fprintf(2,'Press any key to continue...\n'); pause;
end

