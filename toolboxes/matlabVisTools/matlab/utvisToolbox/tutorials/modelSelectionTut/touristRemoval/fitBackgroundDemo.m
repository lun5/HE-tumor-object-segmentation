%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File: fitBackgroundDemo.m
%%% Collects data from 3x3 spatial patchs within an image sequence, all
%%% at the same location.  The goal is to be able to estimate the
%%% background intensity from this gray level data over time.
%%% This code first fits mixture models to this data, selects the most
%%% promising one, and then guesses which component in this mixture
%%% distribution corresponds to the background.  As we shall see,
%%% sometimes this guessing is easy, even when the background is
%%% occluded in the image sequence (by tourists :) slightly more than
%%% half the time.

% Prerequisite:
%  - Go through the modelSelectionTutorial

% Dependencies: iseToolbox and utvisToolbox, as usual.

% Last Modified: Nov 2003
% ADJ

clear
close all;
global matlabVisRoot

if length(matlabVisRoot)==0
  matlabVisRoot = '/h/51/jepson/pub/matlab';
  addpath(matlabVisRoot);
end
addpath([matlabVisRoot '/utvisToolbox/tutorials/modelSelectionTut']);

if ~exist('num2strPad')
  dir0 = pwd;
  cd(matlabVisRoot);
  startup;
  cd(dir0);
end

modelTutPath = [matlabVisRoot '/utvisToolbox/tutorials/modelSelectionTut'];
addpath(modelTutPath);
addpath([modelTutPath '/touristRemoval']);
addpath([modelTutPath '/touristRemoval/']);
addpath([modelTutPath '/touristRemoval/background']);

%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%
TRUE = 1==1;
FALSE = ~TRUE;

USE_MOUSE = FALSE;  %% TRUE in order to mouse in an image region.

maxNoiseSigma = 5;  %% An upper bound for the gray-level noise in the background.

nHist = 50;   % Number of histogram bins to use in plots.

%%%%%%%%%%%% Specify shape of spatial sampling region. %%%%%%%%%%%%%%
% If x0,y0 is a selected pixel (see below), then all pixels
% with (x,y) satisfying:
%   x in x0 + [sampleBounds(1), sampleBounds(3)]
%   y in y0 + [sampleBounds(2), sampleBounds(4)]
% will be selected.
% sampleBounds = [0 0 1 1];   % Specifies a 2x2 neighbourhood.
sampleBounds = [-1 -1 1 1]; % Specifies a 3x3 neighbourhood centered on x0,y0

% Select how much output
quiet = 2; % Show fits with different numbers of inliers
%quiet =  1 % Show only the fit selected by penalized likelihood, and
%quiet =  2 % Plot fit but not penalized likelihood stat and don't pause.
%quiet =  3 % Don't plot anything, don't pause.

% Range of models to be tried.
rangeInliers = 1:7;

%%%%%%%%%%%%%% Derived Constants %%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram quantities
binSize = 256.0/nHist;
binCenter = binSize/2.0;
binLeft = [0:nHist] * binSize;

%%%%%%%%%%%%%% Display Image Sequence %%%%%%%%%%%%%%%%%%%%
% Should be in modelSelectionTut/touristRemoval/background/
h = figure(1); clf;
colormap(gray(256));
for frameNum = 0:2:100
  k = 1+frameNum/2;
  fname = ['cropGrab' num2strPad(frameNum, 4) '.pgm'];
  im(:,:,1+frameNum/2) = pgmRead(fname);
  
  figure(1); image(im(:,:,k));
  pause(0.01);
end


%%%%%%% Show median value of gray levels, at each pixel, over image sequence %%%%%%%%%
medIm = median(im, 3);
figure(1); clf;
subplot(1,2,1);
showIm(medIm); 
title('Median image');
varIm = sum((reshape(im, size(im,1)*size(im,2), size(im,3)) - ...
             repmat(medIm(:), 1, size(im,3))).^2, 2)/size(im,3);
subplot(1,2,2);
showIm(reshape(sqrt(varIm),size(im,1), size(im,2)));
title('Variance image');
% Remnants of the tourist remain...
fprintf(2, 'Press any key to continue...'); pause;


% Initialize the mean image to -1 (values < 0 indicate pixels which have
% not been yet been rocessed.
meanIm = -ones(size(im,1), size(im,2));

%%%%%%%%%%%%%% Select a region of pixel %%%%%%%%%%%%%%%%%%%%%%%
if USE_MOUSE
  fprintf(2, 'Mouse in a SMALL region:\n');
  pnts = ginput(2);
  if size(pnts,1) > 0
    pnts = round(pnts);
  else
    pnts = [20 20; 25 25];
  end
else
  pnts = [22     5;...
          25    29];
end
%%%%%%%%%%%%%%  Loop over selected pixels. %%%%%%%%%%%%%%%%%%%%%%%
for y0 = pnts(1,2): pnts(2,2)
  for x0 = pnts(1,1):pnts(2,1)
    
    %%%%%%%%%%%  Get data from a patch around the selected pixel %%%%%%%%%%%
    nx = size(im,2);
    ny = size(im,1);
    cropBox = [ x0 y0 x0 y0] + sampleBounds;
    cropBox = max([cropBox; ones(1, 4)]);
    cropBox = min([cropBox; [nx ny nx ny]]);
    data = cropSequence(im, cropBox);
    data = data(:);
    nSamp = length(data);

    %% Chat to the user while he/she is waiting...
    fprintf(2, '(x0 y0) = (%d, %d)\n', x0, y0);
    if quiet < 3
      %%%%%%%%%%%  Plot data histogram %%%%%%%%%%%
      figure(2); clf;
      set(gca,'FontSize', 16);
      [dataHist binHist]= histo(data, -binSize, binCenter);

      dataHist = dataHist'/nSamp;
      plotData(binHist, dataHist);
      title(sprintf('Data Samples\n(nSamp = %d)', nSamp));
      pause(0.1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% Fit mixture models %%%%%%%%%%%%%%%%%%%%%%%%
    fracCrossVal = 0.0;
    fitTrainResults = zeros(length(rangeInliers), 1);
    for kTrial = 1:length(rangeInliers)

      nInliers = rangeInliers(kTrial);

      %%%%%%%%%%%%%%%%%% Generate a uniform spaced initial guess.  %%%%%%%%%%%%
      mix = ones(1,nInliers+1)/(nInliers+1);
      if nInliers > 0
        tmp = (256/nInliers);
        mn = ([0:(nInliers-1)]+0.5)*tmp;
        sig = ones(1,nInliers) * (tmp/3);
      else
        mn = [];
        sig = [];
      end
      guessModel = struct('nInliers', nInliers,...
                          'mix', mix,...
                          'mean', mn, ...
                          'sig', sig,...
                          'outLike', 1.0/256);

      %%%%%%%%%%%%%%%%%%%%%%%%%  Do EM fit.  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      emOpt.tolLike = 1.0e-01;   % convergence criterion
      emOpt.maxIts = 50;       % Max Number of iterations for EM

      [fitModel{kTrial}, own, totalLogLike, kIts] = fitMixModel(data, ...
                                                  guessModel, emOpt);

      fitTrainResults(kTrial) = totalLogLike(kIts);
      
      if quiet < 1
        %%%%%%%%%%%%%%%%%%  Plot Fit Results %%%%%%%%%%%%%%%%%%%
        
        %% Plot the mixture model and its components.
        pFitComp = zeros(nHist, fitModel{kTrial}.nInliers+1);
        for k=1:nHist
          [pSum p] = probMix(fitModel{kTrial}, binLeft(k:k+1));
          pFitComp(k, :) = p';
        end
        pFitModel = sum(pFitComp,2);

        figure(3); clf;
        subplot(1,2,1);
        set(gca,'FontSize', 16);
        yMax = max(pFitModel)*1.2;  % for scaling the plots
        plotMix([pFitModel pFitComp], yMax);  
        tmp = axis;
        yMax = tmp(4);
        title(sprintf('Fit Model (nInliers = %d)', nInliers));

        %% Plot the mixture model on the same grid used for data
        %% histograms.
        subplot(1,2,2);
        set(gca,'FontSize', 16);
        plotData(binHist, dataHist, yMax);
        hold on;
        plotMix([pFitModel pFitComp], yMax);  
        tmp = axis;
        yMax = tmp(4);
        title(sprintf('Train Data and Fit Model\n(nSamp = %d)', nSamp));
        hold off;

        pause(0.1);
      end  % of plotting fit results

    end % of loop over nInliers

    %%%%%%%%%%%%%% Model Selection %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Penalized Likelihood
    %%%%% If fracCrossVal > 0, you should also try penalized likelihood
    %%%%% on all the data (i.e. with fracCrossVal = 0.0
    priorTrainModel = zeros(length(rangeInliers),1);
    for kTrial = 1:length(rangeInliers)
      prior = priorMixModel(fitModel{kTrial});
      priorTrainModel(kTrial) = log(prior);
    end
    [tmp kVal] = max(fitTrainResults + priorTrainModel);
    mPenalized = fitModel{kVal};

    if quiet < 2
      %% Plot results of the penalized likelihood over different models.
      plotLike(rangeInliers, fitTrainResults + priorTrainModel, []);
      ylabel('log Prob');
      title('Penalized Likelihood'); 
      fprintf(2,...
              'nSamp %d: nInliers %d (penLike)\n',...
              nSamp, mPenalized.nInliers);
      pause(0.5);
    end

    if quiet < 3
      %%%  Plot the Selected Model only.
      pFitComp = zeros(nHist, mPenalized.nInliers+1);
      for k=1:nHist
        [pSum p] = probMix(mPenalized, binLeft(k:k+1));
        pFitComp(k, :) = p';
      end
      pFitModel = sum(pFitComp,2);

      figure(3); clf;
      subplot(1,2,1);
      set(gca,'FontSize', 16);
      yMax = max(pFitModel)*1.2;  % for scaling the plots
      plotMix([pFitModel pFitComp], yMax);  
      tmp = axis;
      yMax = tmp(4);
      title(sprintf('Selected Model\n(nInliers = %d)', mPenalized.nInliers));

      subplot(1,2,2);
      set(gca,'FontSize', 16);
      plotData(binHist, dataHist, yMax);
      hold on;
      plotMix([pFitModel pFitComp], yMax);  
      tmp = axis;
      yMax = tmp(4);
      title(sprintf('Model with Data\n(nSamp = %d)', nSamp));
      hold off;
    end
    
    %% Display the model selected by penalized likelihood.
    mPenalized

    %%%%%%%%%%%%% GUESS background component from mixture model %%%%%%%%%%%%%%%%%%%
    %% Look for the highest peak with a std dev <= maxNoiseSigma gray levels.  The
    %% height of a Gaussian component at the peak is proportional to
    %% mix_k/sigma_k, where mix_k and sigma_k are the mixture proportions
    %% and the standard deviation for the k-th component.
    %% If there is no component with std-dev <= maxNoiseSigma, then
    %% panic... choose the component with the smallest standard deviation.
    backLayer = 0;
    objMax = 0.0;
    for k = 1:mPenalized.nInliers
      obj = mPenalized.mix(k)/mPenalized.sig(k);
      if mPenalized.sig(k) > maxNoiseSigma
        obj = 0;
      end
      if obj > objMax
        backLayer = k;
        objMax = obj;
      end
    end
    if backLayer == 0
      minSig = inf;
      for k = 1:mPenalized.nInliers
        if mPenalized.sig(k) < minSig
          minSig = mPenalized.sig(k);
          backLayer=k;
        end
      end
    end
    
    %% Chat to user, and record mean of selected component (if any
    %% component was selected.
    if backLayer > 0
      fprintf(2, 'Background Layer: %d mean %f\n', backLayer,...
              mPenalized.mean(backLayer));
      meanIm(y0, x0) =  mPenalized.mean(backLayer);
    else
      fprintf(2, 'Unknonw background layer\n');
    end
    
    %% Pause for the user to view selected model and data histogram.
    if quiet < 1
      fprintf(2, 'Press any key to continue...\n');
      pause;
    elseif quiet < 2
      pause(1);
    end


  end
end

%%%%%%%%%%% Display Image Sequence with Background Estimates inserted %%%%%%%%%%%%%
pb = [min(pnts,[],1) - 0.5; max(pnts,[],1)+0.5];
for reps = 1:5
  h = figure(1); clf;
  colormap(gray(256));
  for k = 1:size(im,3)
    backIm = im(:,:,k);
    backIm(meanIm>=0) = meanIm(meanIm>=0);
    figure(1); image(backIm);
    hold on;
    plot([pb(1,1) pb(2,1) pb(2,1) pb(1,1) pb(1,1)], ...
         [pb(1,2) pb(1,2) pb(2,2) pb(2,2) pb(1,2)], 'r');
    title('Fitted background (inset region)');
    pause(0.1);
  end
  pause(1.0);
end

%%% Conclusion:

%  1. There are many examples of a mixture model being fit to 1D data.
%     This was my primary motivation to writing this demo.  Run it with
%     a small patch and quiet = 0 to see all the models that are fit,
%     along with the one that is selected.

%  2. Guessing the background component is often simple, but
%     sometimes it can be dicy.  This isn't ready for prime time.
%     But this approach does illustrate that  a mixture model
%     representation of a data set is general, and can be useful.

%  3. Including  multiple guesses for the background component and a
%     confidence might be interesting.  That is, the background model
%     could itself be a mixture model (eg. see Stauffer and Grimson).

%  4. Much more information is available in the data to improve our guesses:
%      a) the temporal sequence of the data is informative.  Note that
%         currently the mixture model should produce identical results
%         (up to numerical round-off effects) on any permutation of the
%         data.
%      b) local motion information, or spatial neighbours, could be
%         taken into account in guessing the background component.
%     See the literature on background estimation, beginning with
%     Stauffer and Grimson, ICCV, '99.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End fitBackgroundDemo.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
