function cEnt = crossEntMix1D(mixDist, mixTest)
% cEnt = crossEntMix1D(modelDist, modelTest)
% Estimates crossEnt(mixDist||mixTest) == 
%                integral [- mixDist(x) log_e(mixTest(x))]dx

%% All density functions are plotted as histograms on the
%% grid specified by the following parameters
nHist = 500;   % Number of histogram bins
binSize = 256.0/nHist;
binCenter = binSize/2.0;
binLeft = [0:nHist] * binSize;
sqrt2pi = sqrt(2*pi);

% Compute the probability of p(x|G) integrated over a small bin
pDistComp = zeros(nHist, mixDist.nInliers+1);
for k=1:nHist
    [pSum p] = probMix(mixDist, binLeft(k:k+1));
    pDistComp(k, :) = p';
end
%% this is p(x|G)*dx
pGenModel = sum(pDistComp,2);

%% Generate uniform samples for doing integral.
data = [0.5:1:(nHist-0.5)]*binSize;
data = data';

% Data likelihood at uniform samples, given the test mixture model.
[totalLogLike own dataLike] = dataLogLikeMixModel(data, mixTest);
cEnt = -sum(pGenModel .* log(dataLike));

return;
