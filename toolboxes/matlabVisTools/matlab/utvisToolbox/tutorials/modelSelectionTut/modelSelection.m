function [mPenalized, mCrossVal, mLike, trainData, fitModel, seed] = ...
             modelSelection(genModel, nSamp, fracCrossVal, quiet, seed)
%[mPenalized mCrossVal mLike trainData fitModels seed] = ...
%            modelSelection(genModel, nSamp, fracCrossVal, quiet, seed)
%
% Run EM to fit a 1D gray level data set.  Use penalized likelihood,
% cross validation, and maximum likelihood to determine the number of
% components in the Gaussian mixture model.  Also returns all
% models fit in fitModel struct array.
%
% Input (see mixScript.m for useage)
%  genModel - mixture model struct for the generative model.
%  nSamp    - number of pseudo-random data samples to generate from genModel.
% Optional input.
%  fracCrossVal - in [0,1), the fraction of data to use for cross validation,
%                 default 0.
%  quiet    -  2 silent, 1 plots data and genModel only, 0 plots fits too.
%                 default silent.
%  seed     - random number seed for data set generation
%                 default use current random number state.
if (nargin < 3)
  fracCrossVal = 0;
end
if (nargin < 4)
  quiet = 2;  % 2 silent, 1 plots data and genModel only, 0 plots fits too.
end

if (nargin == 5)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Restart the random number generators rand, and randn
  %% given a seed.  This allows us to run  the fits with
  %% and without cross validation using the same data.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rand('state', seed);
  % We also need to start randn using the autmatically generated seedn:
  seedn = round(rand(1,1) * 1.0e+6);
  randn('state', seedn);
end

%%%%%%%%%%% Plotting parameters. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% All density functions are plotted as histograms on the
%% grid specified by the following parameters
nHist = 50;   % Number of histogram bins
binSize = 256.0/nHist;
binCenter = binSize/2.0;
binLeft = [0:nHist] * binSize;

if quiet < 2
  %% Plot the mixture model on the same grid used for data
  %% histograms.
  pGenComp = zeros(nHist, genModel.nInliers+1);
  for k=1:nHist
    [pSum p] = probMix(genModel, binLeft(k:k+1));
    pGenComp(k, :) = p';
  end
  pGenModel = sum(pGenComp,2);

  % pGenComp now contains the component probabilities times the mixing
  % coefficients for the generative model, evaluated on the discrete
  % set of intervals used for histograming.
  % pGenModel is the probability for the overall model (i.e. the sum of
  % the components), evaluated on the same set of intervals.

  %% Plot the generative model 
  figure(1);
  clf;
  subplot(1,2,1);
  set(gca, 'FontSize', 16);
  yMax = max(pGenModel)*1.2;  % for scaling the plots
  plotMix([pGenModel pGenComp], yMax);  
  tmp = axis;
  yMax = tmp(4);
  title('Generative Model');
end


%%%%%%%%%%%% Data Set: Sample the Generative Model %%%%%%%%%%%%%%
data = sampleMix(genModel, nSamp);

if quiet<2
  figure(1);
  subplot(1,2,2);
  set(gca,'FontSize', 16);
  [dataHist binHist]= histo(data', -binSize, binCenter);
  dataHist = dataHist'/nSamp;
  plotData(binHist, dataHist, yMax);
  title(sprintf('Data Samples\n(nSamp = %d)', nSamp));
  %% The plot on the left shows the probabilities of a data
  %% sample landing in one of the histogram bins, according to
  %% the generative model.  The top curve is the mixture distribution
  %% itself, while the other curves are those for each of the 4 components.
  %% Note there is one outlier (flat) component, and 3 Gaussian mixture
  %% components.

  %% The plot on the right is the result of sampling from the generative
  %% model.  Number of samples is nSamp.
  %% Note the histogram has been normalized to sum to one (so the scales
  %% correspond to the plot of the generative model on the left).
  fprintf(2, 'Press any key to continue...\n');
  pause;
end

%%%%%%%  Set the range for the number of inliers in the fit models
maxInliers = 10;
rangeInliers = [0:maxInliers];

%%%%%%%%  Choose a training set and a cross-validation set
fracCrossVal = min([fracCrossVal 1]);
fracCrossVal = max([fracCrossVal 0]);
nCrossVal = ceil(fracCrossVal * nSamp);
nTrain = max([1, nSamp - nCrossVal]);
trainData = data(1:nTrain);
[trainDataHist binHist]= histo(trainData', -binSize, binCenter);
trainDataHist = trainDataHist'/nTrain;

if (fracCrossVal > 0)
  crossValData = data(nTrain+1:end);
  nCrossVal = length(crossValData);
  if nCrossVal > 0
    crossValResults = zeros(length(rangeInliers),1);
  end
else
  crossValData = [];
  crossValResults = [];
  nCrossVal = 0;
end

%%%%%%%%%%%% Loop over the number of inliers in the guess complexity %%%%%%%%
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

  [fitModel{kTrial}, own, totalLogLike, kIts] = fitMixModel(trainData, ...
                                           guessModel, emOpt);

  fitTrainResults(kTrial) = totalLogLike(kIts);
 
  if fracCrossVal > 0
    crossValResults(kTrial) = dataLogLikeMixModel(crossValData,...
                                    fitModel{kTrial});
  end

  if quiet < 1
    %%%%%%%%%%%%%%%%%%  Plot Fit Results %%%%%%%%%%%%%%%%%%%
    %% Plot the mixture model on the same grid used for data
    %% histograms.
    pFitComp = zeros(nHist, fitModel{kTrial}.nInliers+1);
    for k=1:nHist
      [pSum p] = probMix(fitModel{kTrial}, binLeft(k:k+1));
      pFitComp(k, :) = p';
    end
    pFitModel = sum(pFitComp,2);

    figure(2); clf;
    subplot(1,2,1);
    set(gca,'FontSize', 16);
    yMax = max(pFitModel)*1.2;  % for scaling the plots
    plotMix([pFitModel pFitComp], yMax);  
    tmp = axis;
    yMax = tmp(4);
    title(sprintf('Fit Model (nInliers = %d)', nInliers));

    subplot(1,2,2);
    set(gca,'FontSize', 16);
    plotData(binHist, trainDataHist, yMax);
    hold on;
    plotMix([pFitModel pFitComp], yMax);  
    tmp = axis;
    yMax = tmp(4);
    title(sprintf('Train Data and Fit Model\n(nTrain = %d)', nTrain));
    hold off;

    fprintf(2, 'Press any key to continue...\n');
    pause;
  end  % of plotting fit results

end % of loop over nInliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the data total likelihood for each model complexity (i.e. nInliers).
if (quiet < 2)
  plotLike(rangeInliers, fitTrainResults, crossValResults);
end

%% RED LINE: Note that the log likelihood of the fit model typically
%% increases with increasing number of inlier components in the 
%% fit model.  This is because the parameterized models form nested subsets.
%% That is, any model with n components can be represented exactly as a
%% model with m components, when m > n (by setting the extra
%% mixing proportions to zero).  Therefore the true maximum log likelihood
%% MUST increase with increasing number of inlier components.
%% In practice, however, it may not always strictly increase since
%% a local maximum may be found by the EM code instead of the global
%% maximum.  All these fits were done with one simple style
%% of initial guess, with uniformly spaced gaussian components
%% (see above), so finding local maximum could be an issue here.
%% Smarter methods for initial guess generation, like random sampling,
%% could be considered.

%% BLUE LINE: When cross-validation is used, ie fracCrossVal > 0
%% Some of the increase in the logLikelihood for the
%% models fit to the training data (red line) is due to over-fitting.
%% How can the machine tell if the model is over-fitting?  
%% One approach, called cross-validation, is to split the available
%% data in half.  Half the data is used to form the training set,
%% and the other half forms the cross-validation set.  We can
%% check for over-fitting by computing the logLikelihood of the
%% fit model on the cross-validation set (which was NOT used in
%% training).  The result for the crossValidation
%% data set is the blue line. 

%% So what do these cross-validation results tell us?
%%
%% Which model has the maximum log-likelihood on a
%% cross-validation data set?  Should we select this
%% as the appropriate model?  (Note we cannot select
%% the model with the maximum likelihood, since this
%% will typically be the model with the most components,
%% or nearly so, and it may be massively over-fitted.)

%% Below we consider selecting the model with the maximum cross-validation
%% log likelihood.

%%%%%%%%%%%%%% Model Selection %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Maximum Likelihood 
[tmp kVal] = max(fitTrainResults);
mLike = fitModel{kVal};

%%%%% Cross Validation
if nCrossVal > 0
  [tmp kVal] = max(crossValResults);
  mCrossVal = fitModel{kVal};
else
  mCrossVal = [];
end

if (quiet < 2 & fracCrossVal > 0)
  fprintf(2,...
      'nTrain %d: nInliers %d (maxLike), %d (crossVal)\n',...
      nTrain, mLike.nInliers, mCrossVal.nInliers);
  fprintf(2, 'Press any key to continue...\n');
  pause;
end

%%%%% Penalized Likelihood (on TRAINING data only)
%%%%% If fracCrossVal > 0, you should also try penalized likelihood
%%%%% on all the data (i.e. with fracCrossVal = 0.0
priorTrainModel = zeros(length(rangeInliers),1);
for kTrial = 1:length(rangeInliers)
  prior = priorMixModel(fitModel{kTrial});
  priorTrainModel(kTrial) = log(prior);
end
[tmp kVal] = max(fitTrainResults + priorTrainModel);
mPenalized = fitModel{kVal};

%% Plot results over training data set.
if quiet < 2
  plotLike(rangeInliers, fitTrainResults + priorTrainModel,...
           crossValResults);
  ylabel('log Prob');
  if nCrossVal > 0
    title('Penalized Likelihood and Cross Validation'); 
  else
    title('Penalized Likelihood'); 
  end
  fprintf(2,...
      'nTrain %d: nInliers %d (maxLike), %d (penLike)\n',...
      nTrain, mLike.nInliers, mPenalized.nInliers);
  fprintf(2, 'Press any key to continue...\n');
  pause;
end

return;
