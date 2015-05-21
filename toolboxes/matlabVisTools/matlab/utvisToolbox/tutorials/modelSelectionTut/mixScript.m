%% FILE: mixScript.m
%% Model Selection tutorial.
%% 
%% Use the EM-algorithm to fit a mixture model to synthetic
%% gray-level image data.  Compare three model selection 
%% criteria, namely maximum likelihood, cross validation, and
%% penalized likelihood (approx. MAP).

%% Dependencies:
%%   dataLogLikeMixModel.m modelSelection.m      priorMixModel.m       
%%   fitMixModel.m         plotData.m            probMix.m             
%%   gauss1D.m             plotLike.m            sampleMix.m           
%%   mixScript.m           plotMix.m    
%%         
%% Last Modified: Nov 2003   ADJ

clear
close all;

global matlabVisRoot

% We need to ensure the path is set for the iseToolbox.
if isempty(matlabVisRoot)
  dir = pwd;
  cd ~jepson/pub/matlab   % CHANGE THIS
  startup;
  cd(dir);
end

addpath([matlabVisRoot '/utvisToolbox/tutorials/modelSelectionTut']);

nSamp = 500;  % Total number of data samples to use
fracCrossVal = 0.5;  % Fraction of the data to use for cross validation
                     % For no cross validation, set to 0.

%%%%%%%%%%%  Initialize random number generator %%%%%%%%%%%%%%%%%%%%%%%
% Random number generator seed:
seed = round(sum(1000*clock));
seed0 = seed;

%%%%%%%%%%%% Generative Model %%%%%%%%%%%%%%%%%
%% Define a mixture model with nInliers Gaussian
%% component distributions, and one uniform outlier distribution.
%% Assume the range of the values is [0, 256).
nInliers = 3;
mix = ones(1,nInliers+1)/(nInliers+1);
mn = [64 112 192];
sig = [16 16 16];
genModel = struct('nInliers', nInliers,...
                  'mix', mix,...
                  'mean', mn, ...
                  'sig', sig,...
                  'outLike', 1.0/256);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run several trials of model fitting using different data sets.
nTrain = round((1.0-fracCrossVal) * nSamp);
counts = zeros(11,5);
for kTrial = 1:20
  %% 20 trials MIGHT be enough to see the general trends, but either rerun
  %% mixScript a few times, or use 100 trials instead, to be sure.

  %%%%%%% Fit models using a part of the data, and do model selection %%%%%%%%%
  %% YOU SHOULD OPEN modelSelection.m AND FOLLOW THE COMMENTS THERE.
  if (fracCrossVal > 0.0)
    [mPenalized mCrossVal mLike data dummy seed] = ...
        modelSelection(genModel, nSamp, fracCrossVal, kTrial-1, seed);
    if (nTrain > 0)
      fprintf(2,...
      '%d: nTrain %d: nInliers %d (maxLike), %d (crossVal), %d (penLike)\n',...
       kTrial, nTrain, mLike.nInliers, mCrossVal.nInliers,...
       mPenalized.nInliers);
    end

    %% Keep a count of the number of layers selected by each criteria.
    layers = mLike.nInliers + 1;
    counts(layers, 1) = counts(layers, 1) + 1;
    layers = mCrossVal.nInliers + 1;
    counts(layers, 2) = counts(layers, 2) + 1;
    layers = mPenalized.nInliers + 1;
    counts(layers, 3) = counts(layers, 3) + 1;
  end

  %%%%%%% Fit models using all of the data, and do model selection %%%%%%%%%
  %%%%%%% We cannot do cross validation here.                      %%%%%%%%%
  [mPenalized mCrossVal mLike data dummy seed] = ...
      modelSelection(genModel, nSamp, 0.0, kTrial-1, seed);

  fprintf(2, '   nTrain %d: nInliers %d (maxLike), %d (penLike)\n',...
    nSamp, mLike.nInliers, mPenalized.nInliers);

  layers = mLike.nInliers + 1;
  counts(layers, 4) = counts(layers, 4) + 1;
  layers = mPenalized.nInliers + 1;
  counts(layers, 5) = counts(layers, 5) + 1;
 

  %% Reset random number seed for next iteration
  seed = round(1.0e+6 * rand);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summarize model selection results. %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
set(gca, 'FontSize', 18);
hl=plot(0:10, counts(:,1), 'o-',...
        0:10, counts(:,2), '+-',...
        0:10, counts(:,3), 's-',...
        0:10, counts(:,4), '*-',...
        0:10, counts(:,5), 'd-', 'LineWidth', 2 );
legend(hl, sprintf('MaxLike(%d)', nTrain), sprintf('CrossVal(%d)', nTrain),...
    sprintf('PenLike(%d)', nTrain), sprintf('MaxLike(%d)', nSamp), ...
    sprintf('PenLike(%d)', nSamp));
title('Model Comparison Results');
ylabel('Frequency');
xlabel('nInliers');
% This shows the selected number of inliers given the fitted models to the
% training data and the whole data set. (The number in brackets in the legend
% indicates how many data items were fit.  For the default settings
% above, we used nSamp = 500 and fracCrossVal = 0.5.  So 250 and 500
% data samples were used in the EM-fits.)
% Model selection criteria:
%  MaxLike - Model with the maximum likelihood on the training data
%            (eg. MaxLike(250)), or on the whole data set (eg. MaxLike(500)).
%  CrossVal - Model with the maximum likelihood on the cross validation
%             set.
%  PenLike  - Model with the maximum penalized likelihood, trained on
%             the same data as the cross validation models (eg. PenLike(250)),
%             or trained on the whole data set (eg. PenLike(500)).

%% Conclusion:
%% Maximum likelihood: 
%%   Clearly insufficient in controlling the model
%%   complexity (i.e. nInliers).  The occurence of models with nInliers less
%%   than the maximum allowable value (i.e. 10) indicates cases in which only
%%   local maxima were found by the EM algorithm (or it was stopped
%%   at the maximum number of iterations before converging).
%% Cross validation: 
%%   Manages to limit the model complexity.  Models with more components
%%   are sometimes selected.  Notch at nInliers = 4 might be related to
%%   the particular initial guess.
%% Penalized Likelihood (approximate MAP): 
%%   The method of choice for this problem.
%%   Works well for nSamp about 500.  Get better results when we
%%   use all the data to train the model (i.e. PenLike(500)) instead
%%   of just half the data (i.e. PenLike(250)).

%%%%%%%%%%  Report initial random seed for restarting %%%%%%%%%%%%%%%
fprintf(2, 'Initial random seed: %d\n', seed0);
%%%%%%%%%%  
fprintf(2,'Press any key to continue...');
pause; fprintf(2,'ok\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Entropy and Cross Entropy %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The entropy of a distribution p(x) is 
%%   int_x -p(x) log(p(x)) dx  
%% where int_x denotes the itegral over x.
%% Similarly the cross entropy of q(x) wrt p(x) is
%%   crossEnt(p, q) = int_x -p(x) log(q(x)) dx
%% The cross entropy of q(x) (wrt a fixed p(x)) is minimized 
%% iff q(x) = p(x).

%% Let's check out the cross entropies of fitted models.

%% Fit models using all of the data 
seed = round(sum(1000*clock));
[mPenalized mCrossVal mLike trainData fitModel seed] = ...
      modelSelection(genModel, nSamp, 0.0, kTrial-1, seed);

%% Compute the cross entropy of all the fitted models with the true
%% distribution.  
%% For later use we also compute the log data likelihood and
%% the prior for all these models
clear cEnt layers dataFitLogLike prior;
for k = 1:length(fitModel)
  cEnt(k) = crossEntMix1D(genModel, fitModel{k});
  layers(k) = fitModel{k}.nInliers;
  dataFitLogLike(k) = dataLogLikeMixModel(trainData, fitModel{k});
  prior(k) = priorMixModel(fitModel{k});
end

%% Compute the entropy of the generative model.
entGen = entMix1D(genModel);

%% Plot the cross entropies of the fitted models, and 
%% compare with the entropy of the generative model.
figure(1); clf; 
set(gca, 'FontSize', 18);
h = [];
h(1) = plot(layers, cEnt, 'bo-', 'LineWidth',2);
hold on;
h(2) =plot([0 layers(end)], [entGen entGen], 'r-','LineWidth',2);
legend(h, 'Fit Models', 'Generative Model');
title('Cross Entropy of Fitted Models');
ylabel('Entropy (nats)');
xlabel('nInliers');
%% The minimum entropy fit might occur for nInliers = 3.
%% NOTE: the cross entropies are always the same or larger than the 
%% entropy of the generative model.
fprintf(2,'Press any key to continue...');
pause; fprintf(2,'ok\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Cross Entropy and Expected Data Log Likelihood %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Given a fit model distribution q(x) and the generative distribution
%% p(x) for the data. What is the expected data log likelihood, according
%% to the model distribution q(x), for N samples from p(x)? 
%%   Expected data log Likelihood = N E[log(q(x))]
%% where the expectation E[.] is taken over the data distribution p(x).
%% In terms of the cross entropy of q wrt p, this works out to be:
%%    Expected data log Likelihood = -N crossEnt(p, q). 

%% Compare the expected data log likelihood to the one actually fit
%% on a sample of the data.  .
figure(2); clf; 
set(gca, 'FontSize', 18);
h = [];
h(1) = plot(layers, dataFitLogLike, 'ro-', 'LineWidth',2);
hold on;
h(2) =plot(layers, -length(trainData)*cEnt, 'bo-','LineWidth',2);
legend(h, 'Train Data Log Like', '-N * crossEnt', 0);
title('Log Likelihood of Fitted Models');
ylabel('Log Like');
xlabel('nInliers');
fprintf(2,'Press any key to continue...');
pause; fprintf(2,'ok\n');


%% If we sample sets of data from the genrative model, of length
%% N=nSamp, and compute the data log likelihood according to our 
%% current fitted models, then the mean of these should be the
%% blue curve in the previous figure.
dataLogLike = [];
for t=1:10
  data = sampleMix(genModel, nSamp);
  for k = 1:length(fitModel)
    dataLogLike(k) = dataLogLikeMixModel(data, fitModel{k});
  end
  figure(2); hold on;
  h(3) =plot(layers, dataLogLike, 'go-','LineWidth',1);
end
legend(h, 'Train Data Log Like', '-N * crossEnt', 'Test Data Log Like', 4);
fprintf(2,'Press any key to continue...');
pause; fprintf(2,'ok\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Expected Data Log Likelihood and Model Comparison %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The log unnormalized posterior for a mixture model q(x)
%% is given by 
%%    logUnPost(q) = log q(D) + log(prior(q))
%% From above the expected value of this, over different
%% random data sets of length N, is
%%  E[logUnPost(q)] =  -N crossEnt(p,q) + log(prior(q))
%% For a given fitted distribution q, this is just a straight
%% line in N, the size of the data sets D.

%% We will plot these lines next.  To resolve the lines for the
%% different q's it is better to plot the difference between the expected
%% log un-normalized posterior and N ent(p):
%%  E[logUnPost(q)] - N ent(p) = N (ent(p) - crossEnt(p,q))+log(prior(q))
%% By subtracting out N ent(p) we are removing the common trend
%% in all the data log-likelihoods.
%% For a given value of N, the spacing between the various
%% E[logUnPost(q)] is not affected by removing this trend in N.
col = hsv(length(fitModel));
figure(1); clf; 
set(gca, 'FontSize', 18); hold on;
nData = 0:100:1000;
for k = 1:length(fitModel)
  h(k) = plot(nData, nData*(entGen - cEnt(k)) +log(prior(k)), ...
              'Color', col(k,:), 'LineWidth',2);
  names(k, :) = sprintf('%2d inliners', fitModel{k}.nInliers);
end
set(gca, 'XTick', 0:200:1000);
axis([0 1000 -75 5]); grid on;
title('Bayesian Model Comparison');
ylabel('Log Unnormalized Posterior');
xlabel('nData');
set(gca, 'FontSize', 12);
legh = legend(h, names,-1);
fprintf(2,'Press any key to continue...');
pause; fprintf(2,'ok\n');

%% In this Bayesian Model Comparison plot notice:
%%
%%   1) The log-prior shows up as the y-intercept of the lines.
%%
%%   2) Simpler models (fewer inliers) have larger priors.
%%
%%   3) The lines L_k(N) have negative (or zero) slopes (k==nInlier):
%%         m = -(crossEnt(p,q_k) - ent(p))
%%      where p is the generative model, and q_k one of the fit models.
%%      The term:
%%         D(p || q_k) = crossEnt(p,q_k) - ent(p)
%%      is the Kullback-Leibler (KL) divergence for densities p and q_k.
%%      Since the cross-entropy is always equal to or larger than the
%%      entropy, D(p || q_k) is always non-negative. Also, D(p || q_k) = 0
%%      only when p and q_k are the same density.
%%
%%   4) Models q_k which are poor approximations of the generative 
%%      model p have large KL divergences, and hence large negative
%%      slopes (eg. k == nInlier = 0 or 1).
%%
%%   5) Better approximations q of the distribution p have smaller
%%      KL divergences, and hence flatter slopes (eg. nInlier >=3).
%%      (If you increase nSamp to 5000, and restart at the 
%%      Entropy and CrossEntropy comment above, you will find the
%%      models that are fit for nInlier>=3 more closely approximate p,
%%      and hence the slopes of their lines in this Bayesian Model
%%      Comparison plot are shallower).
%%
%%   6) In Bayesian model comparison we fit various models
%%      q_k(x) to some given sampled data.  The generative density
%%      p(x) for this data is unknown.
%%      If we have set aside a test data set of N samples, then the Bayesian
%%      model comparison plot, for this value of N, shows the relative
%%      sizes of the log unnormalized posterior that can be expected
%%      by evaluating this test data.  The posterior probability for
%%      model k can then be approximated from this plot as
%%         post_k  = exp(L_k(N))/[sum_j L_j(N)]
%%      In particular, the top line in this plot for any N is the
%%      model that is most likely to be chosen using this criterion.
%%      Moreover, the separation between the lines tells us the
%%      ratio of the posteriors,  
%%            post_k/post_j = exp(L_k(N) - L_j(N)).
%%      So a separation L_k(N) - L_j(N) = -5 in this plot signifies
%%      that q_j is exp(-5) = 0.0067 times as probable, given
%%      this data, than q_k.
%%
%%   7) Whether or not we can resolve the 3 inlier components in
%%      this data depends largely on how much data we have, and on
%%      the priors.   
%%      Figure 1 shows that for 600 or more data points, the 3-inlier
%%      model q_3 is the clear cut winner.  For somewhere between 100 and 300
%%      data points, we might expect the 2-inlier model to be selected
%%      as the most probable.  For fewer than 50 data points, the 0-inlier
%%      model narrowly wins.  This result stands so long as the models
%%      q_k(x) are trained on enough data (eg. nSamp>=250).
%%
%%   8) Remember that the individual models q_k(x) were fit to
%%      the sample data (of size nSamp).  Figure 1 shows
%%      the expected log-unPosterior when these models are
%%      evaluated on an independent set of data of size N,
%%      as is done in cross validation.  However, the simple 
%%      penalized likelihood method evaluates the log-unPosterior on
%%      the training data, NOT on a cross-validation set. 
%%      The log data likelihood on the training data, per data sample, is just 
%%          dataFitLogLike(k)/nSamp
%%      We can use this in place of the negative cross entropy.
col = hsv(length(fitModel));
figure(1); hold on;
set(gca, 'FontSize', 18); hold on;
nData = 0:100:1000;
for k = 1:length(fitModel)
  h(k) = plot(nData, ...
              nData*(entGen + dataFitLogLike(k)/nSamp) +log(prior(k)), ...
              'Color', col(k,:), 'LineWidth',1, 'LineStyle', '--');
  names(k, :) = sprintf('%2d inliners', fitModel{k}.nInliers);
end
%%      The resulting dashed lines illustrate the penalized likelihood
%%      method without cross-validation.  Notice that these
%%      dashed lines can slope upwards.  The point on the line
%%      at N = nSamp denotes the penalized likelihood attributed
%%      to each model...the ordering and spacing at this N determines the
%%      posterior prob's assigned to various interpretations,as above.
%%      Notice that, for N around nSamp (and nSamp >= 500) the relative spacing
%%      and ordering of these dashed lines is roughly the same as for the
%%      solid lines (the expected value of the cross validation results).
%%      However, for smaller nSamp (eg nSamp = 100), there can be a
%%      significant difference between the cross entropy and the average
%%      data fit log likelihood (i.e. between the slopes of the solid 
%%      and dashed lines).  This difference is due to two effects: a)
%%      the variation of the data log likelihood around its expected
%%      value; and b) the bias in the data log likelihood on the training
%%      set due to the fitting procedure itself (often the red curve
%%      in Figure 2 is above the blue curve).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Things to do:
%% 1. Try changing nSamp.  
%%    Eg. nSamp = 200. From the plots of the data it is often hard to see
%%    evidence for three inlier components.  Intuitively a model with
%%    two gaussians often appears to be sufficient.  Indeed, the penalized 
%%    likelihood model selects 2 inlier components more often.
%%    Eg. nSamp = 1000 (or more).  Penalized likelihood is almost
%%    perfect in deciding on 3 inlier components.  Cross validation
%%    does not improve significantly from the performance with nSamp=500.
%%
%% 2. Try changing the generative model.
%%    The left two Gaussian components in the current generative model
%%    are harder to resolve, especially for smaller values of nSamp.
%%    If you increase their standard deviations slightly, and/or move
%%    them closer, they will require more data to resolve them.
%%    Alternatively, decreasing their standard deviations and/or moving
%%    them further apart, will make the inference problem easier.

