
%% Dependencies: iseToolbox, utvisToolbox.
%% Data files: eyes.mat nonEyes.mat kevinDwn.pgm

TRUE = (1==1);
FALSE = ~TRUE;
useDisplay = TRUE;
sizeIm = [25 20];

%%%%%%%%%%%  Initialize random number generator %%%%%%%%%%%%%%%%%%%%%%%
% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;

%%%%%%%%%%%%%% Load eye and non-eye images %%%%%%%%%%%%%%%%%%%%%%%%%%
load('eyes');
trainIm0 = [leftEyes rightEyes];
nEyePairs = min([size(leftEyes,2) size(rightEyes,2)]);

load('nonEyes');
nNonEyes = size(nonEyes, 2);

%%%%%%%%%%%%%% Display some of the eye images %%%%%%%%%%%%%%%%%%%%%%%%%%
%% The original face images were scaled and rotated so that, 
%% in the warped image, the left and right eyes have a horizontal separation
%% of 28 pixels.  These warped images were cropped to 20x25 regions
%% centered on each eye.  leftEyes contains the left eye of the person
%% being imaged (on the right side of the picture for you), and
%% similarly for rightEyes.
figure(1); clf;
for k=1:20
  b = 1 + floor(nEyePairs * rand(1,1));
  subplot(1,2,1); 
  showIm(reshape(rightEyes(:,b), sizeIm), [0 255]);
  subplot(1,2,2);
  showIm(reshape(leftEyes(:,b), sizeIm), [0 255]);
   fprintf(2, 'Eye images %d.  ', b);
  % fprintf(2, 'Mean coeff %f\n', meanImU' * [rightEyes(:,b) leftEyes(:,b)]);
  fprintf(2, 'Press any key to continue...\n');
  pause;
end
%% Notice the change in contrast of the various images, along
%% with changes of the person, the lighting, position, pose
%% (eg open/closed eyes), and the occasional reflection from
%% glasses.

%%%%%%%%%%%
%% Rescale training set and identify outliers.
%%%%%%%%%%%
trainIm = trainIm0;
%%
%% Since PCA involves a least squares estimate, it is
%% important to remove outliers from the training data before
%% processing.  As we have seen, least squares methods can
%% be quite sensitive to large outliers.  (The outliers
%% in the current data set are not such a big problem, but
%% we'll go through the exercise of removing them.)
%%
%% A similar issue is the relative scaling of the training 
%% data.  If some images have much larger variance than others,
%% then these images will dominate the others in the computation
%% of the principle vectors.  Images with a range of different
%% contrasts will have the property that the images with
%% higher contrasts will have a larger variance.
%% 
%% Contrast is defined in terms of the amount of variation
%% of the image intensities compared to the average brightness.
%% We estimate the contrast next.
%%
%% To determine the average brightness, we use a constant "DC" image
dcIm = ones(prod(sizeIm),1)/sqrt(prod(sizeIm));
% dcIm is a constant image normalized to unit length.
dcIm' * dcIm

%% To estimate the amount of brightness variation in an image (due
%% to the underlying signal, not to noise) we will use
%% the following mean image (with the constant term removed).
meanIm = sum(trainIm, 2)/size(trainIm,2);
% Remove DC from meanIm
meanIm = meanIm - (dcIm' * meanIm) * dcIm;
% Rescale mean image to be a unit vector.
meanIm = meanIm /sqrt(sum(meanIm .* meanIm));
% Get the std dev of the variation in meanIm, per pixel
meanSig = sqrt(sum(meanIm .* meanIm)/prod(sizeIm))
% Since we have scaled meanIm to be a unit vector, the maximum
% and minimum values are small.
[mn mx] = range2(meanIm)

%% Show the mean image, scaled to be a unit vector.
figure(1); clf;
showIm(reshape(meanIm, sizeIm));

%% Compute components in dcIm and meanIm directions
dcCoeff = dcIm' * trainIm;  % These are the average image brightnesses
meanCoeff = meanIm' * trainIm; 

%% Plot variance per pixel
varCoeff = sum((trainIm - dcIm * dcCoeff) .^ 2, 1)/prod(sizeIm);
figure(1); clf;
plot(meanCoeff, varCoeff, '.b');
figAxis = axis;
hold on;
%% Plot the component of the variance due to the meanIm variation.
plot(0:figAxis(2), ([0:figAxis(2)] * meanSig).^2, 'g');
figAxis(1) = min([figAxis(1) 0]);
axis(figAxis);
ylabel('Gray level Variance (/pixel)');
xlabel('Mean Coeff');
title('Total variation of eye images');
%% Notice that there are a large range of variances in the training
%% images, with those images with larger values of the meanCoeff
%% having significantly larger variances.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Balance the variance in the elements of the training set by
%% projecting out the mean image and adjusting for contrast.

%% Remove DC and meanIm components from training set
trainIm = trainIm - dcIm * (dcIm' * trainIm);
trainIm = trainIm - meanIm * (meanIm' * trainIm);

%% Plot residual variance in training set
varCoeff = sum(trainIm .^ 2, 1)/prod(sizeIm);
figure(1); clf;
plot(meanCoeff, varCoeff, '.b');
ylabel('Variance (/pixel)');
xlabel('Mean Coeff');
title('Component of Variance w/o Mean Image');

%% There is still a contrast dependent trend.
figure(1); clf;
plot(meanCoeff, varCoeff.^0.5, '.');
ylabel('Std. Dev. (/pixel)');
xlabel('Mean Coeff');
title('Component of Std.Dev. w/o Mean Image');

%% Either the spatial structure of the eye images change
%% as the mean component increases, or they stay roughly
%% the same and only the variation of these components
%% increases with the mean component.  Indeed, if the
%% variation of the mean component is primarily due to
%% lighting and imaging effects, then we might assume that
%% the underlying signal is invariant to contrast.
%%
%% Assuming the latter, we attempt to reweight the data set
%% to a constant contrast.
%%
%% Fit a line through the origin to the std dev vs mean component data.
%% (This is biased a bit by outliers...)
A = [meanCoeff; varCoeff.^0.5]';
size(A)
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

%% In order to reweight the data, it is convenient to assume
%% a minimum value for the variance, say due to indendent pixel noise,
%% and add this to the above estimate.
minSig = 4;
plot(xRange(1):xRange(2),...
     sqrt(minSig^2 + (pixelSigRate*[xRange(1):xRange(2)]).^2) , 'r');
%% We will rescale the data so that the variance of points on
%% the red curve will be one.

%% Rescale training set by mean coefficient.
rescalePar = [minSig pixelSigRate];
rescaleCoeff = (minSig^2 + (meanCoeff * pixelSigRate).^2).^0.5;
rescaleIm = trainIm ./ repmat(rescaleCoeff, prod(sizeIm), 1);

%% Plot residual variance in rescaled training set
varCoeff = sum(rescaleIm .^ 2, 1)/prod(sizeIm);
figure(1); clf;
plot(meanCoeff, varCoeff .^ 0.5, '.');
xlabel('Mean Coeff');
ylabel('Std. Dev. (per pixel)');
title('Rescaled Image Std.Dev.'); 
% The majority of the images have been rescaled to have
% a variance of roughly one.  
hold on; ax=axis; plot([ax(1) ax(2)], [1 1], 'r');

%%%%%%%%%%%%% Extreme images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% There are a small number images with large variances and/or negative
%% mean component.  These might be due to outliers in the data set.
%% Let's take a look at them. Given this scatter plot of variances
%% for the rescaled images, we define an extreme image as one with
%% a pixel variance of at least 4^2, or a mean image coefficient
%% less than 20.
trimVar = 4.0 ^ 2.0;
minMeanCoeff = 40;
%% Notice that this mean coefficient corresponds to a mean image
%% with the following range of gray levels
[mn mx] = range2(meanIm * minMeanCoeff)
%% That it, over the whole eye region, it only varies by +- 4 or 5
%% gray levels!  
figure(1); clf; showIm(reshape(meanIm * minMeanCoeff, sizeIm));

%%%%%%%% Trim the extreme images (and view them)
rescaleCoeff = (minSig^2 + (meanCoeff * pixelSigRate).^2).^0.5;
rescaleIm = trainIm ./ repmat(rescaleCoeff, prod(sizeIm), 1);
ok = (sum(rescaleIm .^ 2, 1)/prod(sizeIm) <= trimVar) & ...
     (meanCoeff > minMeanCoeff);
trimIndices = find(~ok);
figure(1); clf;
for b = trimIndices;
  showIm(reshape(trainIm0(:,b), sizeIm), [0 255]);
  fprintf(2, 'Eye images %d.\n', b);
  fprintf(2, 'Mean coeff: %f, Residual std.dev. %f\n', ...
          meanCoeff(b), sqrt(varCoeff(b)));
  fprintf(2, 'Press any key to continue...\n');
  pause;
end
%% We have thrown out this many images
sum(~ok)
%% ... in percent:
sum(~ok) /size(trainIm0,2) * 100
%% In most of these there are easily identifiable problems, such as 
%% the eye is closed, not centered, or at a different scale, or there
%% are reflections from glasses. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Throw out trimmed images before SVD.
trainIm = trainIm0(:, ok);
nTrain = size(trainIm,2);

%% Recompute mean image and rescaling...
%% Compute component in dc direction
dcCoeff = dcIm' * trainIm;
%% Recompute mean image.
meanIm = sum(trainIm-dcIm*dcCoeff, 2)/size(trainIm,2);
% Rescale
meanIm = meanIm /sqrt(sum(meanIm .* meanIm));
figure(1); clf;
showIm(reshape(meanIm, sizeIm));

%% Compute component in meanIm direction
meanCoeff = meanIm' * trainIm;
%% Rescale training set by mean coefficient.
rescalePar = [minSig pixelSigRate];
rescaleCoeff = (minSig^2 + (meanCoeff * pixelSigRate).^2).^0.5;
rescaleIm = (trainIm - dcIm * dcCoeff - meanIm * meanCoeff)./...
              repmat(rescaleCoeff, prod(sizeIm), 1);

%% Plot residual variance in training set
varCoeff = sum(rescaleIm .^ 2, 1)/prod(sizeIm);
figure(1); clf;
scatter(meanCoeff, varCoeff .^ 0.5, '.');
xlabel('Mean Coeff');
ylabel('Std. Dev. (per pixel)');
title('Rescaled Training Image Std.Dev.'); 

%%%%%%%
% Do the SVD on the trimmed rescaled training eyes
%%%%%%%
if (size(rescaleIm,2) > size(rescaleIm,1))
  [v s basisIm] = svd(rescaleIm',0);
else
  [basisIm s v] = svd(rescaleIm,0);
end
sigmaUniEye = diag(s)/sqrt(size(rescaleIm,2));
clear v s
sigma = sigmaUniEye;

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
%% Note the rapid drop off in the singular values.  The first
%% singular direction explains about 20% of the total variance,
%% the second about 8 or 9% and so on... as follows:
dQ([1:5]) * 100

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
%% Note the first five singular directions account for almost half
%% the variance, as follows:
Q([1:5]) * 100

% Plot log(dQ(k)) over a large range of components
% (all but the last two... which give artificially small
% or zero variances).
x = 1:(size(sigmaUniEye,1)-2);
figure(1); clf;
plot(x, log(sigma2(x)/totalVar));
xlabel('Singular value index');
ylabel('Log_e(Fraction of Variance)');
title('Variance Fraction Explained by Subspace of Dimension k');
%% Note the singular values drop off rapidly at first, then
%% decrease roughly exponentially (linearly in this plot) for
%% a large range of k values.  This is typical for a wide
%% variety of natural stimuli.

%%%%%%%%%%%%%%
%%% Plot the first nBasis basis images
%%%%%%%%%%%%%% 
nBasis = 50;
for k=1:nBasis
  figure(1); clf;
  showIm(reshape(basisIm(:,k), sizeIm));
  title(sprintf('Basis Image #%d',k));
  fprintf(1,'Press any key to continue...\n');
  pause;
end
%% Notice how the basis images for small values of k depict
%% larger scales, while the scale (of blobs of one sign, for example)
%% gets finer as k increases.


%%%%%%%%%%%%%%
%%% Orthogonality of basis images.
%%%%%%%%%%%%%% 
%% The constant image, the mean image and the basis images with
%% non-zero singular values form an orthonormal set of 500-dimensional vectors
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

%%%%%%%%%%%%%%
%%% Coefficient Histograms
%%%%%%%%%%%%%%
%% Histograms of coefficients are very roughly Gaussian, mean zero,
%% and standard deviation sigmaUniEye(k) for component k
%% A strong assumption of the detection model is that these
%% distributions can be reasonably described in terms of Gaussians.
nRescale = size(rescaleIm, 2);
figure(1); clf;
superimpose=FALSE;  % superimpose=FALSE; % shows individual histograms
for cDim=1:50
  coeffs = basisIm(:,cDim)' * rescaleIm/sigmaUniEye(cDim);
  figure(2); clf;
  plot(coeffs);

  [hist bin] = histo(coeffs, 64);
  histDens = hist/nRescale/(bin(2)-bin(1));
  % This normalizes the density so that:
  %    sum( histDens ) * (bin(2)-bin(1))
  % is one.
  figure(1);
  if superimpose
     hold on;
  end
  plot(bin, histDens);
  axis([-8 8 0 0.8]);
  hold on;
  dx = 0.1;
  x = [-6:dx:6];
  y = 1.0/sqrt(2.0*pi) * exp(-(x .* x)/2);
  plot(x, 1.0/sqrt(2.0*pi) * exp(-(x .* x)/2), 'g');
  % Here sum(y) * dx is one.
  lam = sqrt(2.0);
  plot(x, lam/2.0 * exp(-abs(x*lam)), 'r');
  xlabel('coeff(k)/sigma(k)');
  ylabel('Probability');
  title(sprintf('Coeff histogram for comp. # %d',cDim));
  hold off;
  pause;
end
%% From the plot the Gaussian approx appears to be
%% a reasonable model, although the data is certainly
%% not clear on this.  For example, a Laplace distribution
%% might be at least as appropriate.

%% Show 2D scatter plots of the eigen image coefficients.
useMean = TRUE;
sigma = sigmaUniEye;
for cDim = 2:10
  figure(1); clf;
  plot(basisIm(:,1)' * rescaleIm, basisIm(:,cDim)' * rescaleIm, 'or'); 
  grid on;
  axis equal;
  title(sprintf('Regular eyes projected onto components 1 and %d', cDim));
  xlabel('Coeff. of 1st comp.');
  ylabel(sprintf('Coeff for comp. # %d', cDim));
  hold on;
  if useMean 
    theta = [-1:0.01:1]*pi;
  else
    theta = [-1/2:0.01:1/2]*pi;
  end
  x = sigma(1)*cos(theta);
  y = sigma(cDim) * sin(theta);
  plot(x,y,'g', 2*x, 2*y, 'b');
  hold off;
  cDim
  fprintf(1, 'Type any key to continue...\n');
  pause;
end

%%%%%%%%%%%%%%%%%%% FILL IN CODE FOR QUESTION 1 HERE %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Detector Characteristics %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FILL IN CODE FOR QUESTION 2 HERE 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Eigen-Eye Detector  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read a test image:
im = pgmRead('kevin.pgm');
clf; showIm(im);

%% FILL IN CODE FOR QUESTION 3 HERE %%%%%%%%%%%%%%%%%%%%

