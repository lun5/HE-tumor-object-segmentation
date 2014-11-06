%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  File: detectEigenEyes.m
%  Matlab script file
%  Date: Oct, 03
%
% Use the PCA model for eye images trained in trainEigenEyes to detect eyes.
% Both out-of-subspace and within-subspace statistics are considered,
% although the out-of-subspace errors appear to be the most significant.

% Dependencies, Toolboxes:
%      iseToolbox/
%      utvisToolbox/file/
% Data file: trainAndTestEyes.mat nonEyes.mat

% Things still to do:
%   Robust fit.
%   SPCA

% Author: ADJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Check Path and Constants  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

%  Check path.  If these aren't on your path, then you will need
%  to add these toolboxes to your path.  See ~jepson/pub/matlab/startup.m
which showIm          % should be iseToolbox\pyrTools\showIm.m

TRUE = (1==1);
FALSE = ~TRUE;


%%%%%%%%  Initialize random number generator %%%%%%%%%%%%%%%%%%%%%%%
% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;

%%%%%%%%%%% Load eye and non-eye images %%%%%%%%%%%%%%%%%%%%%%%%%%
load('trainAndTestEyes');
nTest = size(testIm,2);
nTrain = size(trainIm,2);

load('nonEyes');
nNonEyes = size(nonEyes, 2);

% Constant image
dcIm = ones(prod(sizeIm),1)/sqrt(prod(sizeIm));

% Do detection on the testIm (not the trainIm, used to train the PCA model)
eyeIm = testIm;
%eyeIm = trainIm;  % Check overfitting.

%%%%%%%%%%%%%%%%%
% Project out dc and mean image

dcNonEyeCoeff = dcIm' * nonEyes;
meanNonEyeCoeff = meanIm' * nonEyes;
testNonEyeIm = nonEyes - dcIm * dcNonEyeCoeff - meanIm * meanNonEyeCoeff;

dcEyeCoeff = dcIm' * eyeIm;  
meanEyeCoeff = meanIm' * eyeIm; 
testEyeIm = eyeIm - dcIm * dcEyeCoeff - meanIm * meanEyeCoeff;


%%%%%%%%%% Residual Variance Plots and ROC curves %%%%%%%%%%%%%%%%%%%%%%%%%
trial_nBasis = [0 1 2 4 8 20 50 100 200];
thetaThres = zeros(size(trial_nBasis));
meanThres = 20;
USE_VAR_IM = TRUE;
USE_SQR_ERROR = TRUE;

col = hsv(length(trial_nBasis));
figure(2); clf;
figure(3); clf;
h1 = [];
h2 = [];
legendNames = [];
for kTrial = 1:length(trial_nBasis)
  nBasis = trial_nBasis(kTrial);
  
  if USE_VAR_IM
    varIm=basisIm(:,(nBasis+1):end).^2 * sigmaUniEye((nBasis+1):end).^2;
  else
    varIm = ones(prod(sizeIm), 1);
  end

  % Compute residual images and in-space statistics
  if nBasis > 0  % Project out first nBasis directions
    coeffs = basisIm(:,1:nBasis)' * testEyeIm;
    residualEyeIm = testEyeIm - basisIm(:,1:nBasis) * coeffs;

    coeffs = basisIm(:,1:nBasis)' * testNonEyeIm;
    residualNonEyeIm = testNonEyeIm - basisIm(:,1:nBasis) * coeffs;
  
  else 
    residualEyeIm = testEyeIm;
    varEyeIns = zeros(1, size(testEyeIm, 2));
    residualNonEyeIm = testNonEyeIm;
    varNonEyeIns = zeros(1, size(testNonEyeIm, 2));
  end
  
  % COmpute out of subspace statistic.
  if USE_SQR_ERROR
    varCoeff = repmat(varIm, 1, size(residualEyeIm,2));
    varEyeOos = sum((residualEyeIm.^2)./varCoeff, 1)/prod(sizeIm);
  else 
    varCoeff = repmat(varIm.^0.5, 1, size(residualEyeIm,2));
    varEyeOos = sum(abs(residualEyeIm)./varCoeff, 1)/prod(sizeIm);
  end

  if USE_SQR_ERROR
    varCoeff = repmat(varIm, 1, size(residualNonEyeIm,2));
    varNonEyeOos = sum((residualNonEyeIm.^2)./varCoeff, 1)/prod(sizeIm);
  else
    varCoeff = repmat(varIm.^0.5, 1, size(residualNonEyeIm,2));
    varNonEyeOos = sum(abs(residualNonEyeIm)./varCoeff, 1)/prod(sizeIm);
  end

  % Plot residual variance in non eye
  figure(1); clf;
  hold on;
  if USE_SQR_ERROR
    plot(meanNonEyeCoeff, varNonEyeOos .^ 0.5, 'r.');
    plot(meanEyeCoeff, varEyeOos.^0.5, 'g.');
  else
    plot(meanNonEyeCoeff, varNonEyeOos, 'r.');
    plot(meanEyeCoeff, varEyeOos, 'g.');
  end
  title(sprintf('Out of subspace error, nBasis = %d', nBasis));
  xlabel('Mean Coeff');
  ylabel('Out of subspace error (graylevels/pixel)');

  % Plot ROC curve
  thetaSamp = 0.0:0.001:(pi/2);
  legendNames = [legendNames; sprintf('nBasis =%3d', nBasis)];
  eyeStat = varEyeOos;
  nonEyeStat = varNonEyeOos;
  if USE_SQR_ERROR
    eyeStat = eyeStat .^ 0.5;
    nonEyeStat = nonEyeStat .^ 0.5;
  end
  
  theta = atan2(eyeStat, meanEyeCoeff);
  theta(meanEyeCoeff < meanThres) = inf;
  cummEye = sum(repmat(theta(:)',length(thetaSamp),1) < ...
                repmat(thetaSamp(:),1,length(theta)), 2)/length(theta);
  theta = atan2(nonEyeStat, meanNonEyeCoeff);
  theta(meanNonEyeCoeff < meanThres) = inf;
  cummNonEye = sum(repmat(theta(:)',length(thetaSamp),1) < ...
                repmat(thetaSamp(:),1,length(theta)), 2)/length(theta);
 
  figure(3);  hold on;
  h1(kTrial) = plot(thetaSamp, cummEye, 'Color', col(kTrial, :));
  plot(thetaSamp, 1-cummNonEye, '--', 'Color', col(kTrial, :));
  axis([0 0.4 0.8 1]);
  grid on;
  legend off;
  legend(h1, legendNames, 4);
  title('True Detection (solid) and True Rejection (dashed)');
  xlabel('Theta (rad)');
  ylabel('Rate');
  
  figure(2); hold on;
  h2(kTrial) = plot(cummNonEye, cummEye, 'Color', col(kTrial, :));
  axis([0 0.2 0.8 1]);
  grid on;
  legend off;
  legend(h2, legendNames, 4);
  title('ROC');
  xlabel('False Positive Rate');
  ylabel('True Detection Rate');

  % Plot threshold line for roughly 5% false positives 
  % on out of subspace error.
  kT = find(cummNonEye <= 0.05);
  if length(kT) > 0
    kT = kT(end);
    thetaThres(kTrial) = thetaSamp(kT);
    theta = thetaThres(kTrial);
    x = cos(theta);
    y = sin(theta);
    figure(1); hold on;
    plot([0 1000*x], [0 1000*y], 'k');
    fprintf(2,...
    'nBasis: %3d, True Detection Rate %5.3f, False Detection Rate %5.3f\n',...
            nBasis, cummEye(kT), cummNonEye(kT));
  end
  
  if kTrial == 1
    fprintf(2, 'Separate the three figures\n');
  end
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end  
% The optimal detector (for roughly equal false positive and false
% negative rates) is for:
%    nBasis:  50, True Detection Rate 0.943, False Detection Rate 0.049


% Residual Variance Plots and ROC curves with both out-of subspace
% and within-subspace statistics
figure(3); close 3;
trial_nBasis = [20 50]; % Or use them all: [0 1 2 4 8 20 50 100 200];
thetaThres = zeros(length(trial_nBasis), 2);
USE_VAR_IM = TRUE;
USE_SQR_ERROR = TRUE;

for kTrial = 1:length(trial_nBasis)
  nBasis = trial_nBasis(kTrial);
  
  if USE_VAR_IM
    varIm=basisIm(:,(nBasis+1):end).^2 * sigmaUniEye((nBasis+1):end).^2;
  else
    varIm = ones(prod(sizeIm), 1);
  end

  % Compute residual images and in-space statistics
  if nBasis > 0  % Project out first nBasis directions
    coeffs = basisIm(:,1:nBasis)' * testEyeIm;
    residualEyeIm = testEyeIm - basisIm(:,1:nBasis) * coeffs;
    if USE_SQR_ERROR
      varEyeIns = ...
        sum((coeffs./repmat(sigmaUniEye(1:nBasis),1,size(coeffs,2))).^2,1) ...
        /nBasis;
    else
      varEyeIns = ...
        sum(abs(coeffs./repmat(sigmaUniEye(1:nBasis),1,size(coeffs,2))),1) ...
        /nBasis;
    end

    coeffs = basisIm(:,1:nBasis)' * testNonEyeIm;
    residualNonEyeIm = testNonEyeIm - basisIm(:,1:nBasis) * coeffs;
    if USE_SQR_ERROR
      varNonEyeIns = ...
        sum((coeffs./repmat(sigmaUniEye(1:nBasis),1,size(coeffs,2))).^2,1) ...
        /nBasis;
    else
      varNonEyeIns = ...
        sum(abs(coeffs./repmat(sigmaUniEye(1:nBasis),1,size(coeffs,2))),1) ...
        /nBasis;
    end
  else 
    residualEyeIm = testEyeIm;
    varEyeIns = zeros(1, size(testEyeIm, 2));
    residualNonEyeIm = testNonEyeIm;
    varNonEyeIns = zeros(1, size(testNonEyeIm, 2));
  end
  
  % COmpute out of subspace statistic.
  if USE_SQR_ERROR
    varCoeff = repmat(varIm, 1, size(residualEyeIm,2));
    varEyeOos = sum((residualEyeIm.^2)./varCoeff, 1)/prod(sizeIm);
  else 
    varCoeff = repmat(varIm.^0.5, 1, size(residualEyeIm,2));
    varEyeOos = sum(abs(residualEyeIm)./varCoeff, 1)/prod(sizeIm);
  end

  if USE_SQR_ERROR
    varCoeff = repmat(varIm, 1, size(residualNonEyeIm,2));
    varNonEyeOos = sum((residualNonEyeIm.^2)./varCoeff, 1)/prod(sizeIm);
  else
    varCoeff = repmat(varIm.^0.5, 1, size(residualNonEyeIm,2));
    varNonEyeOos = sum(abs(residualNonEyeIm)./varCoeff, 1)/prod(sizeIm);
  end

  % Plot residual variance in non eye
  figure(1); clf;
  if nBasis > 0
    subplot(1,2,1); hold on;
    if USE_SQR_ERROR
      plot(meanNonEyeCoeff, varNonEyeIns .^ 0.5, 'r.');
      plot(meanEyeCoeff, varEyeIns .^ 0.5, 'g.');
    else
      plot(meanNonEyeCoeff, varNonEyeIns, 'r.');
      plot(meanEyeCoeff, varEyeIns, 'g.');
    end
    title(sprintf('In-space error, nBasis = %d', nBasis));
    xlabel('Mean Coeff');
    ylabel('In space error');
    
  end
  subplot(1,2,2); hold on;
  if USE_SQR_ERROR
    plot(meanNonEyeCoeff, varNonEyeOos .^ 0.5, 'r.');
    plot(meanEyeCoeff, varEyeOos.^0.5, 'g.');
  else
    plot(meanNonEyeCoeff, varNonEyeOos, 'r.');
    plot(meanEyeCoeff, varEyeOos, 'g.');
  end
  title(sprintf('Out of subspace error, nBasis = %d', nBasis));
  xlabel('Mean Coeff');
  ylabel('Out of subspace error (graylevels/pixel)');

  thetaSamp = 0.0:0.001:(pi/2);
  figure(2); clf;

  nA = 6;
  col = hsv(nA);
  figure(2); clf;
  h1 = [];
  h2 = [];
  legendNames = [];
  if nBasis==0
    kaRange = 1;
  else
    kaRange = 1:nA;
  end
  for ka = kaRange
    a = (ka-1) * 1/(nA-1);
    legendNames = [legendNames; sprintf('a = %4.2f', a)];
    eyeStat = (a * varEyeIns + (1-a)*varEyeOos);
    nonEyeStat = (a * varNonEyeIns + (1-a)*varNonEyeOos);
    if USE_SQR_ERROR
      eyeStat = eyeStat .^ 0.5;
      nonEyeStat = nonEyeStat .^ 0.5;
    end
    
    % Plot ROC curve
    
    theta = atan2(eyeStat, meanEyeCoeff);
    theta(meanEyeCoeff < meanThres) = inf;
    cummEye = sum(repmat(theta(:)',length(thetaSamp),1) < ...
                  repmat(thetaSamp(:),1,length(theta)), 2)/length(theta);
    theta = atan2(nonEyeStat, meanNonEyeCoeff);
    theta(meanNonEyeCoeff < meanThres) = inf;
    cummNonEye = sum(repmat(theta(:)',length(thetaSamp),1) < ...
                     repmat(thetaSamp(:),1,length(theta)), 2)/length(theta);
    
    figure(2); 
    subplot(1,2,1);
    hold on;
    h1(ka) = plot(thetaSamp, cummEye, 'Color', col(ka, :));
    plot(thetaSamp, 1-cummNonEye, '--', 'Color', col(ka, :));
    axis([0 0.4 0.8 1]);
    grid on;
    legend off;
    legend(h1, legendNames, 4);
    title('True Detection/Rejection Rates');
    xlabel('Theta (rad)');
    ylabel('Rate');

    
    figure(2); subplot(1,2,2); hold on;
    h2(ka) = plot(cummNonEye, cummEye, 'Color', col(ka, :));
    axis([0 0.2 0.8 1]);
    grid on;
    legend off;
    legend(h2, legendNames, 4);
    title('ROC');
    xlabel('False Positive Rate');
    ylabel('True Detection Rate');

  end
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end 
% The optimal detector uses a approximately 0.  Using just the
% out-of-subpace error corresponds to a = 0.  Therefore the
% within-subspace statistic is not adding much to the detection
% performance (at least in this style of detection).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Check Out False Responses %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nBasis = 20;
meanThres = 20;
USE_MOUSE = TRUE;
USE_VAR_IM = TRUE;
USE_SQR_ERROR = TRUE;
USE_MAX_OOS = FALSE;  % When true, use a maximum threshold on oos error.
maxStat = 40;


% Compute residual images and in-space statistics
if nBasis > 0  % Project out first nBasis directions
  coeffs = basisIm(:,1:nBasis)' * testEyeIm;
  residualEyeIm = testEyeIm - basisIm(:,1:nBasis) * coeffs;

  coeffs = basisIm(:,1:nBasis)' * testNonEyeIm;
  residualNonEyeIm = testNonEyeIm - basisIm(:,1:nBasis) * coeffs;
  
else 
  residualEyeIm = testEyeIm;
  varEyeIns = zeros(1, size(testEyeIm, 2));
  residualNonEyeIm = testNonEyeIm;
  varNonEyeIns = zeros(1, size(testNonEyeIm, 2));
end

% COmpute out of subspace statistic.
if USE_VAR_IM
  varIm=basisIm(:,(nBasis+1):end).^2 * sigmaUniEye((nBasis+1):end).^2;
else
  varIm = ones(prod(sizeIm), 1);
end

if USE_SQR_ERROR
  varCoeff = repmat(varIm, 1, size(residualEyeIm,2));
  varEyeOos = sum((residualEyeIm.^2)./varCoeff, 1)/prod(sizeIm);
else 
  varCoeff = repmat(varIm.^0.5, 1, size(residualEyeIm,2));
  varEyeOos = sum(abs(residualEyeIm)./varCoeff, 1)/prod(sizeIm);
end

if USE_SQR_ERROR
  varCoeff = repmat(varIm, 1, size(residualNonEyeIm,2));
  varNonEyeOos = sum((residualNonEyeIm.^2)./varCoeff, 1)/prod(sizeIm);
else
  varCoeff = repmat(varIm.^0.5, 1, size(residualNonEyeIm,2));
  varNonEyeOos = sum(abs(residualNonEyeIm)./varCoeff, 1)/prod(sizeIm);
end

eyeStat = varEyeOos;
nonEyeStat = varNonEyeOos;
if USE_SQR_ERROR
  eyeStat = eyeStat .^ 0.5;
  nonEyeStat = nonEyeStat .^ 0.5;
end

% Plot residual variance in non eye
figure(1); clf;
hold on;
plot(meanNonEyeCoeff, nonEyeStat, 'r.');
plot(meanEyeCoeff, eyeStat, 'g.');
title(sprintf('Out of subspace error, nBasis = %d', nBasis));
xlabel('Mean Coeff');
ylabel('Out of subspace error (graylevels/pixel)');

kTrial = find(trial_nBasis == nBasis);
if USE_MOUSE || size(kTrial) == 0
  
  fprintf(2,'Mouse in endpoint of radial line for threshold in Fig. 1\n');
  pt = ginput(1);
  theta0 = atan2(pt(1,2), pt(1,1));

  if USE_MAX_OOS
    maxStat = pt(2);
  end
  
else
  theta0 = thetaThres(kTrial);
end

% Plot threshold line.
figure(1); hold on;
x = cos(theta0);
y = sin(theta0);
figure(1); hold on;
plot([0 1000*x], [0 1000*y], 'k');
if USE_MAX_OOS
  % Start at s0*[x y] = [x0 maxStat], i.e. s0 = maxStat/y
  plot([(maxStat/y)*x 1000*x], [maxStat maxStat],'k');
end

theta = atan2(eyeStat, meanEyeCoeff);
theta(meanEyeCoeff < meanThres) = inf;
if USE_MAX_OOS
  theta(eyeStat > maxStat) = inf;
end
detEye = theta < theta0;
theta = atan2(nonEyeStat, meanNonEyeCoeff);
theta(meanNonEyeCoeff < meanThres) = inf;
if USE_MAX_OOS
  theta(nonEyeStat > maxStat) = inf;
end
detNonEye = theta < theta0;

fprintf(2,'True detection rate: %5.3f (%d errors out of %d)\n', ...
       sum(detEye)/length(detEye), length(detEye)-sum(detEye), length(detEye));
fprintf(2,'False detection rate: %5.3f  (%d errors out of %d)\n', ...
       sum(detNonEye)/length(detNonEye), sum(detNonEye), length(detNonEye));

FN = ~detEye;
FP = detNonEye;

% Show false positives
for kIm = find(FP(:)')
  
  im = nonEyes(:,kIm);
  figure(2); clf;
  subplot(1,2,1);
  showIm(reshape(im, sizeIm), [0 255]);
  title(sprintf('FP %d/%d, nonEye %d', sum(FP(1:kIm)), sum(FP), kIm));
  
  d = dcIm' * im;
  m = meanIm' * im;
  if nBasis > 0
    c = basisIm(:, 1:nBasis)' * im;
  else
    c = [];
  end
  recon = d * dcIm + m * meanIm;
  if nBasis > 0
    recon = recon +  basisIm(:, 1:nBasis) * c;
  end
  subplot(1,2,2);
  showIm(reshape(recon, sizeIm), [0 255]);
  title(sprintf('Recon nBasis=%d', nBasis));
  
  
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end
% This detector is still accepting a lot of images which look nothing
% like eyes!!!!!

% Show false negatives
for kIm = find(FN(:)')
  
  im = eyeIm(:,kIm);
  figure(2); clf;
  subplot(1,2,1);
  showIm(reshape(im, sizeIm), [0 255]);
  title(sprintf('FN %d/%d, Eye %d', sum(FN(1:kIm)), sum(FN), kIm));
  
  d = dcIm' * im;
  m = meanIm' * im;
  if nBasis > 0
    c = basisIm(:, 1:nBasis)' * im;
  else
    c = [];
  end
  recon = d * dcIm + m * meanIm;
  if nBasis > 0
    recon = recon +  basisIm(:, 1:nBasis) * c;
  end
  subplot(1,2,2);
  showIm(reshape(recon, sizeIm), [0 255]);
  title(sprintf('Recon nBasis=%d', nBasis));
  
  
  fprintf(2, 'Press any key to continue...');
  pause; fprintf(2, 'ok\n');
end
% This detector is rejecting a lot of eye images which look
% like reasonable eyes!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Further Experiments %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Try using the same images for training and testing.  See
% the comment above, "Check overfitting" and set eyeIm = trainIm.
% Then rerun the subsequent code to compute the ROC curves.  Notice
% that for nBasis = 100 or more, this reports better detection performance
% than when the test set was used.  This suggest we are using
% information specific to the training set to do the detection.

% Try rerunning the detection only but with a maximum
% allowable out of subspace statistic.  (i.e. set 
% USE_MAX_OOS = TRUE; USE_MOUSE=TRUE;
% and set the threshold in the previous section.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% End: Detect Eigen-Eye %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
