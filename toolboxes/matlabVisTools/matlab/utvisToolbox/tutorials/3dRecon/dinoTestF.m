% File: dinoTestF
% A4 2003 handout code
% Estimate F matrix from synthetic corresponding points.
%
% ADJ, Nov. 03

clear
close all
FALSE = 1 == 0;
TRUE = ~FALSE;

sigmaNoise = 1.0;  % Std dev (in pixels) of image point locations.

global matlabVisRoot % Directory containing ise- and utvis-toolboxes.

% We need to ensure the path is set for the iseToolbox.
if isempty(matlabVisRoot)
  dir = pwd;
  cd /h/51/jepson/pub/matlab   % CHANGE THIS to your startup directory
  startup;
  cd(dir);
end

reconRoot = [matlabVisRoot '/utvisToolbox/tutorials/3dRecon'];  
addpath([reconRoot '/utils']);


% Random number generator seed:
ranSeed = round(sum(1000*clock));
rand('state', ranSeed);
ranSeed0 = ranSeed;
% We also need to start randn. Use a seedn generated from seed:
ranSeedn = round(rand(1,1) * 1.0e+6);
randn('state', ranSeedn);

nTrial = 10;  % Number of ransac trials to use


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up cameras
% The cameras are automatically rotated by projectDino to fixate
% on the mean of the 3D points.  We do not always want to allow
% the cameras to move in this fashion (eg. when we change sclZ).
% So we will compute the rotations for the left and right cameras
% once and for all, and then use these.
f = 100; % focal length
dLeft = [-50, 0, -150]';  % Center of projection for left camera
dRight = [50, 0, -150]';  % Center of projection for right camera
% Compute camera rotations to fixate on Dino's center.
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, [], 1.0);
Rleft = MextLeft(:, 1:3);
[pRight polys MintRight MextRight] = projectDino(f, dRight, [], 1.0);
Rright = MextRight(:, 1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate data...
sclZ = 1.0;
% Dino left image data
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, Rleft, sclZ);

% Dino right image data
[pRight polys MintRight MextRight] = projectDino(f, dRight, Rright, sclZ);

% Add image noise
pLeft(1:2, :) = pLeft(1:2, :) + sigmaNoise * randn(2,size(pLeft,2));
pRight(1:2, :) = pRight(1:2, :) + sigmaNoise * randn(2,size(pRight,2));

% Show left and right images
hFig = figure(1); clf; 
plot(pLeft(1,:), pLeft(2,:), '.b');
axis xy; axis equal;
xlabel('X'); ylabel('Y');
axis([-150 150 -100 100]);
title('Left image of Dino');
pause(0.1);

hFig = figure(2); clf; 
plot(pRight(1,:), pRight(2,:), '.b');
axis xy; axis equal;
axis([-150 150 -100 100]);
xlabel('X'); ylabel('Y');
title('Right image of Dino');
pause(0.1);


% Build correspondence data
clear imPts;
imPts = cat(3, pLeft, pRight);
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
hFig = figure(1);
clf; hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,1), imPts(2,:,1), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
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
hFig = figure(2);
clf; hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,2), imPts(2,:,2), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Right Image');
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
