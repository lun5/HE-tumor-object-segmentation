%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  File: trainEigenEyes.m
%  Matlab script file
%  Date: Oct, 03
%
% Train a PCA model for eye images.  The model of a given
% image is in the form
%        im = a * dcIm + b * meanIm + sum_k U_k c_k 
% where dcIm is a constant image, meanIm is a mean eye
% image and U_k, k = 1 to nBasis are PCA directions.  The
% form of this model allows for detecting eyes over a range of
% contrasts.

% Dependencies, Toolboxes:
%      iseToolbox/
%      utvisToolbox/file/
% Data file: eyes.mat

% Things still to do:
%   Train multiple PCA patches (multiple means)
%   Effect of blurring image data.
%   SPCA

% Author: ADJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Check Path and Constants  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

%  Check path.  If these aren't on your path, then you will need
%  to add these toolboxes to your path.  See ~jepson/pub/matlab/startup.m
which showIm          % should be iseToolbox\pyrTools\showIm.m
which histo           % should be iseToolbox\pyrTools\MEX\histo.*

TRUE = 1 == 1;
FALSE = ~TRUE;

sizeIm = [25 20];

%%%%%%%%  Initialize random number generator %%%%%%%%%%%%%%%%%%%%%%%
% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;

%%%%%%%%%%% Load eye and non-eye images %%%%%%%%%%%%%%%%%%%%%%%%%%
load('eyes');
eyeIm0 = [leftEyes rightEyes];
nEyePairs = min([size(leftEyes,2) size(rightEyes,2)]);
nEyes0 = size(eyeIm0,2);

%  Select half the eye images for training, reserve the other half for
%  testing.
idx = randperm(nEyes0);
nEyes = floor(nEyes0/2);
eyeIm = eyeIm0(:, idx(1:nEyes));
testIm = eyeIm0(:, idx((nEyes+1):end));
nTest = size(testIm,2);

%%%%%%%%%%% Display some of the eye images %%%%%%%%%%%%%%%%%%%%%%%%%%
% The original face images were scaled and rotated so that, 
% in the warped image, the left and right eyes have a horizontal separation
% of 28 pixels.  These warped images were cropped to 20x25 regions
% centered on each eye.  leftEyes contains the left eye of the person
% being imaged (on the right side of the picture for you), and
% similarly for rightEyes.
figure(1); clf;
for k=1:20
  b = 1 + floor(nEyePairs * rand(1,1));
  subplot(1,2,1); 
  showIm(reshape(rightEyes(:,b), sizeIm), [0 255]);
  title(sprintf('Right Eye %d', b));
  subplot(1,2,2);
  showIm(reshape(leftEyes(:,b), sizeIm), [0 255]);
  title(sprintf('Left Eye %d', b));
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end
% Notice the change in contrast of the various images, along
% with changes of the person, the lighting, position, pose
% (eg open/closed eyes), and the occasional reflection from
% glasses.

%%%%%%%%%%%
% Rescale training set and identify outliers.
%%%%%%%%%%%
%
% Since PCA involves a least squares estimate, it is
% important to remove outliers from the training data before
% processing.  As we have seen, least squares methods can
% be quite sensitive to large outliers.  (The outliers
% in the current data set are not such a big problem, but
% we'll go through the exercise of removing them.)
%
% A similar issue is the relative scaling of the training 
% data.  If some images have much larger variance than others,
% then these images will dominate the others in the computation
% of the principle vectors.  Images with a range of different
% contrasts will have the property that the images with
% higher contrasts will have a larger variance.
% 
% Contrast is defined in terms of the amount of variation
% of the image intensities compared to the average brightness.
% We estimate the contrast next.
%
% To determine the average brightness, we use a constant "DC" image
dcIm = ones(prod(sizeIm),1)/sqrt(prod(sizeIm));
% dcIm is a constant image normalized to unit length.
dcIm' * dcIm

% To estimate the amount of brightness variation in an image (due
% to the underlying signal, not to noise) we will use
% the following mean image (with the constant term removed).
meanIm = sum(eyeIm, 2)/size(eyeIm,2);
% Remove DC from meanIm
meanIm = meanIm - (dcIm' * meanIm) * dcIm;
% Rescale mean image to be a unit vector.
meanIm = meanIm /sqrt(sum(meanIm .* meanIm));
% Get the std dev of the variation in meanIm, per pixel
meanSig = sqrt(sum(meanIm .* meanIm)/prod(sizeIm))
% Since we have scaled meanIm to be a unit vector, the maximum
% and minimum values are small.
[mn mx] = range2(meanIm)

% Show the mean image, scaled to be a unit vector.
figure(1); clf;
showIm(reshape(meanIm, sizeIm));
title('Mean Eye Image');
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

% Compute components in dcIm and meanIm directions
dcCoeff = dcIm' * eyeIm;  % These are the average image brightnesses
meanCoeff = meanIm' * eyeIm; 

% Plot variance per pixel
varCoeff = sum((eyeIm - dcIm * dcCoeff) .^ 2, 1)/prod(sizeIm);
figure(1); clf;
plot(meanCoeff, varCoeff, '.');
figAxis = axis;
hold on;
% Plot the component of the variance due to the meanIm variation.
plot(0:figAxis(2), ((0:figAxis(2)) * meanSig).^2, 'g');
figAxis(1) = min([figAxis(1) 0]);
axis(figAxis);
ylabel('Gray level Variance (/pixel)');
xlabel('Mean Coeff');
title('Total variation of eye images');
% Notice that there are a large range of variances in the training
% images, with those images with larger values of the meanCoeff
% having significantly larger variances.
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Balance the variance in the elements of the training set by
% projecting out the mean image and adjusting for contrast.

% Remove DC and meanIm components from training set
trainIm = eyeIm - dcIm * (dcIm' * eyeIm);
trainIm = trainIm - meanIm * (meanIm' * eyeIm);


% Plot residual variance in training set
varCoeff = sum(trainIm .^ 2, 1)/prod(sizeIm);
figure(1); clf;
plot(meanCoeff, varCoeff, '.');
ylabel('Variance (/pixel)');
xlabel('Mean Coeff');
title('Component of Variance w/o Mean Image');
% There is still a contrast dependent trend.
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

% Either the spatial structure of the eye images change
% as the mean component increases, or they stay roughly
% the same and only the variation of these components
% increases with the mean component.  Indeed, if the
% variation of the mean component is primarily due to
% lighting and imaging effects, then we might assume that
% the underlying signal is invariant to contrast.
%
% Assuming the latter, we attempt to reweight the data set
% to a constant contrast.
%
% Fit a line through the origin to the std dev vs mean component data.
% (This is biased a bit by outliers...)
figure(1); clf;
plot(meanCoeff, varCoeff.^0.5, '.');
ylabel('Std. Dev. (/pixel)');
xlabel('Mean Coeff');
title('Component of Std.Dev. w/o Mean Image');

A = [meanCoeff; varCoeff.^0.5]';
[u s v] = svd(A, 0);
pixelSigRate = v(2,1)/v(1,1);
figAxis = axis;
figAxis(1) = min([figAxis(1) 0]);
axis(figAxis);
xRange = figAxis(1:2);
ylabel('Std. Dev./pixel');
xlabel('Mean Coeff');
title('Contrast Dependent Image Variation');
% Here is the fitted line...
hold on; plot(xRange, pixelSigRate * xRange , 'g');
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

% In order to reweight the data, it is convenient to assume
% a minimum value for the variance, say due to independent pixel noise,
% and add this to the above estimate.
minSig = 8;
plot(xRange(1):xRange(2),...
     sqrt(minSig^2 + (pixelSigRate*(xRange(1):xRange(2))).^2) , 'r');
% We will rescale the data so that the variance of points on
% the red curve will be one.

% Rescale training set by mean coefficient.
rescalePar = [minSig pixelSigRate];
rescaleCoeff = (minSig^2 + (meanCoeff * pixelSigRate).^2).^0.5;
rescaleIm = trainIm ./ repmat(rescaleCoeff, prod(sizeIm), 1);

% Plot residual variance in rescaled training set
varCoeff = sum(rescaleIm .^ 2, 1)/prod(sizeIm);
figure(1); clf;
plot(meanCoeff, varCoeff .^ 0.5, '.');
xlabel('Mean Coeff');
ylabel('Std. Dev. (per pixel)');
title('Rescaled Image Std.Dev.'); 
% The majority of the images have been rescaled to have
% a variance of roughly one.  
hold on; ax=axis; plot([ax(1) ax(2)], [1 1], 'r');
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

%%%%%%%%%% Extreme images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are a small number images with large variances and/or negative
% mean component.  These might be due to outliers in the data set.
% Let's take a look at them. Given this scatter plot of variances
% for the rescaled images, we define an extreme image as one with
% a pixel variance of at least 5^2, or a mean image coefficient
% less than 20.
trimVar = 5.0 ^ 2.0;
minMeanCoeff = 20;
% Notice that this mean coefficient corresponds to a mean image
% with the following range of gray levels
[mn mx] = range2(meanIm * minMeanCoeff)
% That it, over the whole eye region, it only varies by +-2 (roughly)
% gray levels!  
figure(1); clf; showIm(reshape(meanIm * minMeanCoeff, sizeIm));
title('Mean Image at Minimum Amplitude');
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

%%%%% Trim the extreme images (and view them)
rescaleCoeff = (minSig^2 + (meanCoeff * pixelSigRate).^2).^0.5;
rescaleIm = trainIm ./ repmat(rescaleCoeff, prod(sizeIm), 1);
ok = (sum(rescaleIm .^ 2, 1)/prod(sizeIm) <= trimVar) & ...
     (meanCoeff > minMeanCoeff);
trimIndices = find(~ok);
figure(1); clf;
for b = trimIndices;
  showIm(reshape(trainIm(:,b), sizeIm));
  title(sprintf('Extreme Eye %d', b));
  fprintf(2, 'Mean coeff: %f, Residual std.dev. %f\n', ...
          meanCoeff(b), sqrt(varCoeff(b)));
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end
% We have thrown out this many images
sum(~ok)
% ... as a fraction of the total number of eye images:
sum(~ok) /size(trainIm,2)
% In most of these there are easily identifiable problems, such as 
% the eye is closed, not centered, or at a different scale, or there
% are reflections from glasses. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Throw out trimmed images before SVD.
trainIm = eyeIm(:, ok);
nTrain = size(trainIm,2);

% Recompute mean image and rescaling...
% Compute component in dc direction
dcCoeff = dcIm' * trainIm;
% Recompute mean image.
meanIm = sum(trainIm-dcIm*dcCoeff, 2)/size(trainIm,2);
% Rescale
meanIm = meanIm /sqrt(sum(meanIm .* meanIm));
figure(1); clf;
showIm(reshape(meanIm, sizeIm));
title(sprintf('Mean Image after discarding Extreme Images'));
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

% Compute component in meanIm direction
meanCoeff = meanIm' * trainIm;
% Rescale training set by mean coefficient.
rescalePar = [minSig pixelSigRate];
rescaleCoeff = (minSig^2 + (meanCoeff * pixelSigRate).^2).^0.5;
rescaleIm = (trainIm - dcIm * dcCoeff - meanIm * meanCoeff)./...
              repmat(rescaleCoeff, prod(sizeIm), 1);

% Plot residual variance in training set
varCoeff = sum(rescaleIm .^ 2, 1)/prod(sizeIm);
figure(1); clf;
plot(meanCoeff, varCoeff .^ 0.5, '.');
xlabel('Mean Coeff');
ylabel('Std. Dev. (per pixel)');
title('Rescaled Training Image Std.Dev.'); 
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

%%%%%%%
% Do the SVD on the trimmed rescaled training eyes
%%%%%%%
fprintf(2, 'Calling svd on a %d x %d matrix...', size(rescaleIm));
if (size(rescaleIm,2) > size(rescaleIm,1))
  [v s basisIm] = svd(rescaleIm',0);
else
  [basisIm s v] = svd(rescaleIm,0);
end
fprintf(2, 'done\n');
sigmaUniEye = diag(s)/sqrt(size(rescaleIm,2));
clear v s
sigma = sigmaUniEye;

%%%%%%%%
% Save the results
%%%%%%%%
if FALSE  % Don't overwrite by default.
  save 'trainAndTestEyes' sizeIm meanIm sigmaUniEye basisIm ...
      rescalePar trainIm testIm;
end


%%%%%%%
% Plot the first 20 singular values: dQ plot
%%%%%%%
% Compute total scaled variance of training set.
sigma2 = sigmaUniEye .* sigmaUniEye;
totalVar = sum(sigma2);
dQ = sigma2/totalVar;
Q = cumsum(dQ);
figure(1); clf;
selectInd = 1:20;
plot(selectInd, dQ(selectInd));
hold on;
plot(selectInd, dQ(selectInd), 'bo');
hold off;
xlabel('Singular value index, k');
ylabel('Fraction of Variance, dQ');
title('dQ(k): Variance Fraction Explained by one s.v.');
% Note the rapid drop off in the singular values.  The first
% singular direction explains about 20% of the total variance,
% the second about 8 or 9% and so on... as follows:
dQ(1:5) * 100
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

% Plot Q(k) = variance fraction explained by first k components
figure(1); clf;
selectInd = 1:20;
plot(selectInd, Q(selectInd));
hold on;
plot(selectInd, Q(selectInd), 'bo');
hold off;
xlabel('Singular value index');
ylabel('Fraction of Variance');
title('Variance Fraction Explained by Subspace');
% Note the first five singular directions account for almost half
% the variance, as follows:
Q(1:5) * 100
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

% Plot log(dQ(k)) over a large range of components
% (all but the last two... which give artificially small
% or zero variances).
x = 1:(size(sigmaUniEye,1)-2);
figure(1); clf;
plot(x, log10(sigma2(x)/totalVar));
xlabel('Singular value index');
ylabel('Log_{10}(Fraction of Variance)');
title('Variance Fraction Explained by Subspace of Dimension k');
% Note the singular values drop off rapidly at first, then
% decrease roughly exponentially (linearly in this plot) for
% a large range of k values.  This is typical for a wide
% variety of natural stimuli.
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

%%%%%%%%%%%%%%
% Plot the first nBasis basis images
%%%%%%%%%%% 
nBasis = 50;
for k=1:nBasis
  figure(1); clf;
  showIm(reshape(basisIm(:,k), sizeIm));
  title(sprintf('Basis Image #%d',k));
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end
% Notice how the basis images for small values of k depict
% larger scales, while the scale (of blobs of one sign, for example)
% gets finer as k increases.


%%%%%%%%%%%%%%
% Orthogonality of basis images.
%%%%%%%%%%% 
% The constant image, the mean image and the basis images with
% non-zero singular values form an orthonormal set of 500-dimensional vectors
% Build the DC image, a constant, unit length image ...
dcIm = ones(size(meanIm))/sqrt(prod(sizeIm));
% By the construction in trainEigenFloat, dcIm and meanIm are
% orthonormal, and therefore the following should be close to
% the 2x2 identity matrix:
[dcIm  meanIm]' * [dcIm meanIm]
% Similarly the basis images are orthogonal to the DC image and
% the mean image.  The following should be nearly zero:
max(max(abs(basisIm(:,sigmaUniEye>1.0e-10)' * [dcIm meanIm])))
% Finally the basis images are orthonormal, so the following
% should be essentially zero:
max(max(abs(basisIm' * basisIm - eye(size(basisIm,2)))))
fprintf(2, 'Press any key to continue...');
pause; fprintf(2, 'ok\n');

%%%%%%%%%%%%%%
% Coefficient Histograms
%%%%%%%%%%%%%%
% Histograms of coefficients are very roughly Gaussian, mean zero,
% and standard deviation sigmaUniEye(k) for component k
% A strong assumption of the detection model is that these
% distributions can be reasonably described in terms of Gaussians.
nRescale = size(rescaleIm, 2);
figure(1); clf;
superimpose=TRUE;  % superimpose=FALSE; % shows individual histograms
USE_LOG = FALSE;
for cDim=1:50
  coeffs = basisIm(:,cDim)' * rescaleIm/sigmaUniEye(cDim);
  figure(2); clf;
  plot(coeffs);

  [hist bin] = histo(coeffs, 64);
  histDens = hist/nRescale/(bin(2)-bin(1));
  minHistDens = 0.1/nRescale/(bin(2)-bin(1));
  % This normalizes the density so that:
  %    sum( histDens ) * (bin(2)-bin(1))
  % is one.
  figure(1);
  if superimpose
     hold on;
  end
  if ~USE_LOG
    plot(bin, histDens);
    hold on;
    dx = 0.1;
    x = -6:dx:6;
    y = 1.0/sqrt(2.0*pi) * exp(-(x .* x)/2);
    plot(x, 1.0/sqrt(2.0*pi) * exp(-(x .* x)/2), 'g');
    % Here sum(y) * dx is one.
    lam = sqrt(2.0);
    plot(x, lam/2.0 * exp(-abs(x*lam)), 'r');
    axis([-8 8 0 0.8]);
    ylabel('Probability');
  else
    idx = histDens > 0;
    logHistDens = log10(minHistDens) * ones(size(histDens));
    logHistDens(idx) = log10(histDens(idx));
    plot(bin, logHistDens);
    hold on;
    dx = 0.1;
    x = -6:dx:6;
    y = 1.0/sqrt(2.0*pi) * exp(-(x .* x)/2);
    plot(x, log10(1.0/sqrt(2.0*pi) * exp(-(x .* x)/2)), 'g');
    % Here sum(y) * dx is one.
    lam = sqrt(2.0);
    plot(x, log10(lam/2.0 * exp(-abs(x*lam))), 'r');
    axis([-8 8 -4 1]);
    ylabel('log10 Probability');
  end
  xlabel('coeff(k)/sigma(k)');
  title(sprintf('Coeff histogram for comp. # %d',cDim));
  hold off;
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end
% From the plot the tails of the distribution appear to
% be longer than the Gaussian model (green curve).  For example,
% a Laplace distribution (red curve) might be at least as
% appropriate.  Plotting log frequency instead of frequency
% gives a better idea of the tails of the distributions: 
%   Repaste above code using USE_LOG = TRUE;


% Show 2D scatter plots of the eigen image coefficients.
useMean = TRUE;
sigma = sigmaUniEye;
for cDim = [2:10 20 40 80 160 320 450]
  figure(1); clf;
  plot(basisIm(:,1)' * rescaleIm, basisIm(:,cDim)' * rescaleIm, 'or'); 
  axis equal;
  title(sprintf('Regular eyes projected onto components 1 and %d', cDim));
  xlabel('Coeff. of 1st comp.');
  ylabel(sprintf('Coeff for comp. # %d', cDim));
  hold on;
  if useMean 
    theta = (-1:0.01:1)*pi;
  else
    theta = (-1/2:0.01:1/2)*pi;
  end
  x = sigma(1)*cos(theta);
  y = sigma(cDim) * sin(theta);
  plot(x,y,'g', 2*x, 2*y, 'b');
  grid on;
  hold off;
  cDim
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end

% Show reconstructions
idx = randperm(size(testIm,2));
nBasisTry = [5 20 50];
for t = 1:10
  im = testIm(:,idx(t));
  figure(1); clf;
  subplot(2,2,1); showIm(reshape(im, sizeIm), [0 255]);
  title(sprintf('Image %d', idx(t)));
  for j = 1:3
    nBasis = nBasisTry(j);
    dcCoeff = dcIm' * im;
    meanCoeff = meanIm' * im;
    coeffs = basisIm(:,1:nBasis)' * im;
    modelIm = dcCoeff * dcIm + meanCoeff * meanIm + ...
              basisIm(:, 1:nBasis) * coeffs;
    figure(1); subplot(2,2,1+j);
    showIm(reshape(modelIm, sizeIm), [0 255]);
    title(sprintf('Recon, nBasis=%d', nBasis));
  end
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end



% To get an idea of range of possible images that the PCA model
% represents, we show some random images from the PCA model.  That
% is, we chose nBasis Gaussian coefficients with the std devs given
% by sigmaUniEye.

meanCoeff = 350;
minSig = rescalePar(1);
pixelSigRate = rescalePar(2);
rescaleCoeff = (minSig^2 + (meanCoeff * pixelSigRate).^2).^0.5;
nBasisTry = [5 10 20 50 100 200];
for t = 1:50
  % Generate enough random coefficients for all reconstructions
  coeffs = randn(nBasisTry(end), 1);
  coeffs = coeffs .* sigmaUniEye(1:nBasisTry(end));
  for j = 1:6
    % Show the random image using just the first nBasis coeffs.
    nBasis = nBasisTry(j);
    meanCoeff = meanIm' * im;
    rescaleIm = basisIm(:, 1:nBasis) * coeffs(1:nBasis);
    modelIm = meanCoeff * meanIm + rescaleCoeff * rescaleIm;
    [mn mx] = range2(modelIm);
    % Choose dc response to center the range of the image on greylevel 128.
    dcCoeff = (128 - (mx + mn)/2)/dcIm(1);
    modelIm = dcCoeff * dcIm + modelIm;
    figure(1); subplot(2,3,j);
    showIm(reshape(modelIm, sizeIm), [0 255]);
    title(sprintf('Rand nBasis=%d', nBasis));
  end
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% End: Train Eigen-Eyes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
