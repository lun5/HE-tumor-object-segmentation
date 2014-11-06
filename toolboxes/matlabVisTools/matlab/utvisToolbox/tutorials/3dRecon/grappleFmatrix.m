% File: grappleFmatrix
% A4 2003 handout code
% Uses RANSAC to estimate F matrix from corresponding points.
%
% ADJ, Nov. 03

clear
close all
FALSE = 1 == 0;
TRUE = ~FALSE;
global matlabVisRoot

% We need to ensure the path is set for the iseToolbox.
if isempty(matlabVisRoot)
  dir = pwd;
  cd /h/51/jepson/pub/matlab   % CHANGE THIS
  startup;
  cd(dir);
end

reconRoot = [matlabVisRoot '/utvisToolbox/tutorials/3dRecon']; 
addpath([reconRoot '/data/wadham']);
addpath([reconRoot '/utils']);


% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

nTrial = 10;  % Number of ransac trials to try.

% Wadham left image: use  wadham/001.jpg
imPath = 'data/wadham/'; fnameLeft = '001'; 
im = imread([imPath fnameLeft],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imLeft = imDwn;

% Read correspondence data
load data/wadham/corrPnts5
% Wadham right image data/wadham002-5.jpg use for corrPnts2-5 respectively
fnameRight = '005';
im = imread([imPath fnameRight],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imRight = imDwn;

clear imPts;
imPts = cat(3, im_pos1', im_pos2');
nPts = size(imPts,2);
if size(imPts,1) == 2
  imPts = cat(1, imPts, ones(1,nPts,2));
end

% RANSAC for F
seeds = {};
sigma = 2.0; rho = 2;
for kTrial = 1: nTrial
  % Test out F matrix on a random sample of 8 points
  idTest = randperm(nPts);
  nTest = min(8, nPts);
  idTest = idTest(1:nTest);

  % Solve for F matrix on the random sample
  [F Sa Sf] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),1);
  
  % Compute perpendicular error of all points to epipolar lines
  perpErrL = zeros(1,nPts);
  for k = 1:nPts
    lk = imPts(:,k,2)' * F';
    perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
  end
  
  % Detect inliers
  idInlier = abs(perpErrL) < rho*sigma;
  
  % Count inliers
  nInlier = sum(idInlier);
  if nInlier > 20
    % Store sets of sampled points with at least 20 inliers
    seed.id = idTest;
    seed.idInlier = idInlier;
    seed.nInlier = nInlier;
    seed.F = F;
    
    kSeed = length(seeds)+1
    seeds{kSeed} = seed;
  end
end 
% Done RANSAC trials

% Extract best solution
nInliers = zeros(1, length(seeds));
for ks = 1:length(seeds)
  nInliers(ks) = seeds{ks}.nInlier;
end 
[nM ks] = max(nInliers);
nInliers(ks)

%  Refine estimate of F using all inliers.
F = seeds{ks}.F;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;
sum(idInlier)
% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1: 10
  % Fit F using all current inliers
  [F Sa Sf] = linEstF(imPts(:,idInlier,1), imPts(:,idInlier,2),1);
  
  % Compute perpendicular error to epipolar lines
  perpErrL = zeros(1,nPts);
  for k = 1:nPts
    lk = imPts(:,k,2)' * F';
    perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
  end
  idInlier = abs(perpErrL) < rho*sigma;
  nInlier = sum(idInlier)
  
  % If we have the same set of inliers as the previous iteration then stop.
  if all(idInlier == idInlierOld)
    break;
  end
  idInlierOld = idInlier;
end
  
%%%%%%%%%%% Plot results
nTest = 64;  % Number of epipolar lines to plot
nCol = 16;   % Number of different colours to use.
col = hsv(nCol);  % Colour map.

% Random sample the lines to plot
idLines = randperm(nPts);  
idLines = idLines(1:nTest);

% Show left image
SUPERIMPOSE = TRUE;
hFig = figure(2);
clf;
if SUPERIMPOSE
  image(imLeft);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imDwn), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,1), imPts(2,:,1), '.b');
set(gca,'YDir', 'reverse');
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Left Image');
for kl = 1:length(idLines)
  % Plot interest point location corresponding to epipolar line as a "o"
  % in the same colour as the epipolar line.
  k = idLines(kl);
  plot(imPts(1,k,1), imPts(2,k,1), 'o', 'Color', col(mod(k,nCol)+1,:));
  % Plot epipolar line.
  lk = imPts(:,k,2)' * F';
  epk = cropLineInBox(lk(1:2), lk(3), cropBox); 
  set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
end


% Show right image
SUPERIMPOSE = TRUE;
hFig = figure(3);
clf;
if SUPERIMPOSE
  image(imRight);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imDwn), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,2), imPts(2,:,2), '.b');
set(gca,'YDir', 'reverse');
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Left Image');
perpErrR = [];
for kl = 1:length(idLines)
  % Plot interest point location corresponding to epipolar line as a "o"
  % in the same colour as the epipolar line.
  k = idLines(kl);
  plot(imPts(1,k,2), imPts(2,k,2), 'o', 'Color', col(mod(k,nCol)+1,:));
  % Plot epipolar line.
  lk = imPts(:,k,1)' * F;
  epk = cropLineInBox(lk(1:2), lk(3), cropBox); 
  set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
  perpErrR = [perpErrR (lk * imPts(:,k,2)/norm(lk(1:2)))];
end

% Compute perpendicular distance to epipolar lines in left and right images.
perpErrL = [];
for k = 1:nPts
  lk = imPts(:,k,2)' * F';
  perpErrL = [perpErrL (lk * imPts(:,k,1))/norm(lk(1:2))];
end
perpErrR = [];
for k = 1:nPts
  lk = imPts(:,k,1)' * F;
  perpErrR = [perpErrR (lk * imPts(:,k,2)/norm(lk(1:2)))];
end

% Plot a histogram of the perpendicular distances
err = [perpErrL perpErrR];
err = min(err, 10);
err = max(err, -10);
[n b] = histo(err, 64);
figure(4); clf;
plot(b,n);
title('Distance to epipolar line');
xlabel('Error in pixels');
ylabel('Frequency');

% Count inliers
idL = abs(perpErrL)< rho*sigma;
idR = abs(perpErrR) < rho*sigma;
idInlier = idL & idR;
sum(idInlier)
sum(idInlier)/nPts

% save 'Fcorr' F Sa Sf idInlier nInliers
