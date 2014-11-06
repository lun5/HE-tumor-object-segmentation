%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script file: projectiveMassageDino.m 
%
% Demonstrate 3D affine and Euclidean reconstructions from corresponding
% points in perspective projection.
%
% The perspective factorization algorithm is due to Mahamud, Hebert,
% Omori, and Ponce, CVPR, 2001.  I replaced the generalized eigenvalue
% problem with a few optimization steps along the gradient. 
%
% Note, the self calibration follows the algorithm of Polleyfeys, Koch, and van Gool,
% IJCV, 1998, except we search for a nearby rank 3 Q matrix.
%
% TTD:
%  - The self-calibration method is quite sensitive to rescaling.
%    Some thoughts on how to do this are warranted.
%  - We should be solving a CONSTRAINED (and normalized) least squares
%    problem for Q, where the constraint is that Q has rank 3.  See literature?
%  - Turning the self-calibration off is now a problem.  Scaling issue?
%
% ADJ, Nov. 03. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear
global TRUE FALSE;
global matlabVisRoot;

TRUE = 1==1;
FALSE = ~TRUE;

if isempty(matlabVisRoot)
  dir = pwd;
  cd '/h/51/jepson/pub/matlab';  % CHANGE THIS
  %cd 'C:/work/Matlab'
  startup;
  cd(dir);
end
reconRoot = [matlabVisRoot '/utvisToolbox/tutorials'];

addpath([reconRoot '/3dRecon/utils']);
addpath([reconRoot '/3dRecon/projective']);

% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

% Control Flags
DEBUG = FALSE;         % TRUE => Extra output
HOMOG_DEMO = TRUE;     % Show effect of simple homography on 3D data points
NUM_RESCALE = TRUE;    % Use Hartley's rescaling of image data for stability
SELF_CALIBRATE = TRUE;  % Use only self-calibration constraints in metric
                        % estimation. (FALSE => broken, for some reason)
UNDO_RESCALE = FALSE & NUM_RESCALE; % TRUE => Use original scaling of image data 
                                    % during metric upgrade. (Not good.)

% Parameters
sigma = 1.0;     % Std Dev of noise (in pixels) in point locations
tol = 1.0e-9;    % Convergence tolerance for z-estimation in projective fit
nIm = 10;        % Number of data images to use (try 3-10).
                 % nIm = 2 is enough for projective reconstruction, but not
                 % for metric reconstruction algorithm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read  Dino data set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[v f] = getHalfDino;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Place Dino in a fixed 3D pose in world coordinate frame and display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set canonical 3D pose
R0 = [1 0 0; 0 0 1; 0 -1 0];        % Rotate and reflect dino (concave away).
mn0 = [0 0 0]';                   % Place Dino's mean on Z axis.
P0 = R0 * v' + repmat(mn0, 1, size(v,1));
if size(P0,1)==3
  P0 = [P0; ones(1,size(P0,2))];
end
fprintf(2,'Depth statistics...');
imStats(P0(3,:));
nPts = size(P0,2);

% Surface plot of data.
figure(3); clf;
for k = 1:length(f)
  vf = P0(:, f{k});
  patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
end
set(gca,'YDir', 'reverse', 'FontSize', 14);
axis vis3d; axis square; axis equal;
title('Dino Model');
xlabel('X'); ylabel('Y'), zlabel('Z');
fprintf(2,'Rotate this figure.\n');
fprintf(2,'Press any key to continue...');
pause; fprintf(2,'ok\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create some camera locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Columns of dMat give X,Y,Z locations of camera in world coordinate frame.
dMat = [ 0    0    0   50   -70  -50    10   80     2.5    10
         0   50 -100  -20    10    0   -50  -20    10     -10
      -150 -140 -145 -160  -145  -155 -150 -160  -140    -145];
% Focal lengths for each of these 10 cameras.
fMat = [100 125 150 125 125 120 120 120 120 120]';
% Camera rotation will be chosen to fixate Dino's center.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grab nIm Scaled - Orthographic images  (nIm <= length(fMat) == 10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imPts = [];
zPts = [];
M_GT = [];
for k = 1:nIm
  
  d = dMat(:,k);  % Camera center
  f1 = fMat(k);   % Focal length
  
  % Choose rotation to fixate mn0, say
  % That is solve:  R * (mn0 - d) = [0 0 1]';
  R = eye(3);
  R(:,3) = (mn0 - d)/norm(mn0 - d);
  R(:,2) = R(:,2) - (R(:,3)' * R(:,2)) * R(:,3);
  R(:,2) = R(:,2)/norm(R(:,2));
  R(:,1) = R(:,1) - R(:,2:3) * (R(:,2:3)' * R(:,1));
  R(:,1) = R(:,1)/norm(R(:,1));
  R = R';
  if DEBUG
    fprintf(2,'Sanity check:\n   Should be [0 0 Center_Depth]: %f %f %f\n', ...
            (R * (mn0 - d)));
  end

  % Build camera matrix K and image formation matrix M.
  K = diag([f1, f1, 1]);
  M = K * [R -R*d];
  if size(M_GT,1)>0
    M_GT = [M_GT; M];
  else
    M_GT = M;
  end

  % Compute the projected image locations
  p = M * P0;
  
  % Concatenate ground truth depth values over multiple frames.
  if size(zPts,1)>0
    zPts = cat(3, zPts, p(3,:));
  else
    zPts = p(3,:);
  end

  % Concatenate noisy image points over multiple frames.
  p = p ./ repmat(p(3,:),3,1);
  p(1:2,:) = p(1:2,:) + sigma * randn(2,nPts);  % Add noise
  if size(imPts,1)>0
    imPts = cat(3, imPts, p(1:2,:));
  else
    imPts = p(1:2,:);
  end

  % Show image
  figure(10); clf;
  h = plot(p(1,:), p(2,:), '.b');
  set(gca,'YDir', 'reverse', 'FontSize', 14);
  axis([-150 150 -150 150]);
  title(sprintf('Image %d',k));
  fprintf(2,'Press any key to continue...');
  pause;
  fprintf(2,'ok\n');
  
end  % Forming images.


% Reorder the matrices zPts and imPts in the form nIm by nPts
if size(zPts,1) == 1 && size(zPts,2) == nPts && size(zPts,3) == nIm
  zPts = permute(zPts, [3 2 1]);
  zPts = reshape(zPts, nIm, nPts);
end
if size(imPts,1) == 2
  imPts = cat(1, imPts, ones(1, nPts,nIm));
end
if size(imPts,2) == nPts && size(imPts,3) == nIm
  imPts = permute(imPts,[1, 3, 2]);
end
% Set homogeneous 3rd coord in imPts
imPts(3,:,:) = 1;
imPtsSave = imPts;  % Save the image points in case we mess up...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END of data generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rescale image data (attempting to obtain more numerical stability).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Knum = repmat(eye(3), [1 1 nIm]);
if NUM_RESCALE
  %% Rescale for numerical stability
  mn = sum(imPts(1:2,:,:),3)/nPts;
  %mns = reshape(mn, [2 nIm 1]);
  var = sum(sum((imPts(1:2,:,:)-repmat(mn, [1  1 nPts])).^2,3)/nPts, 1);
  % Scale image points so that sum of variances of x and y = 2.
  scl = sqrt(2./var(:));
  % Sanity: varScl =  var .* reshape(scl.^2, [1 nIm]); % Should be all 2's
  
  % Scale so x and y variance is roughly 1, translate so image mean (x,y) is zero.
  Knum(1:2,3,:) = -reshape(mn, [2, 1, nIm]);
  Knum(1:2,:,:) = Knum(1:2,:,:).*repmat(reshape(scl, [1 1 nIm]), [2, 3,1]);
  for kIm = 1:nIm
    imPts(:,kIm, :) = reshape(Knum(:,:,kIm),3,3) * reshape(imPts(:,kIm, :),3,nPts);
  end
  % Sanity check
  % sum(imPts(1:2,:,:),3)/nPts  % Should be [0 0]'
  % sum(sum(imPts(1:2,:,:).^2,3)/nPts,1) % Should be 2.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do Projective Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess for projective depths
z = ones( nIm, nPts);
z = reshape(z, [1 nIm nPts]);

% Do projective factorization iterations
for its=1:50 
  
  % Data matrix
  D = reshape(imPts .* repmat(z, [3, 1, 1]), 3 * nIm, nPts);
  W = sqrt(sum(D .* D, 1)).^-1;
  D = D .* repmat(W, 3*nIm, 1);

  % Factor data matrix
  if size(D,2)<=size(D,1)
    [U S V] = svd(D,0); 
  else
    [V S U] = svd(D',0);
  end
  S = diag(S);
  fprintf(2, 'it %d, sv(5) %e\n', its, S(5));
  if its == 1
    figure(1); clf; 
    hand1 = plot(S, '-*g', 'MarkerSize', 14, 'LineWidth', 2);
    set(gca,'FontSize', 14);
    title('Data Matrix Singular Values'); 
    xlabel('Singular Value Index');
    ylabel('Sigular Value');
    Sstart = S;
  else
    figure(1); clf; 
    hand1=plot(Sstart, '-*g', 'MarkerSize', 14, 'LineWidth', 2);
    hold on;
    hand2=plot(S, '-*r', 'MarkerSize', 14, 'LineWidth', 2);
    set(gca,'FontSize', 14);
    title('Singular Values of Data Matrix');
    xlabel('Singular Value Index');
    ylabel('Singular Value');
    legend([hand1 hand2], 'Iteration 1', sprintf('Iteration %d', its));
  end
  pause(0.1);

  % Extract estimates for camera matrices, Mfac, and projective points, Pfac
  Mfac = U(:,1:4);
  Sfac = S(1:4);
  Pfac = V(:,1:4)';
  
  % Given Mfac, Sfac and Pfac, re-estimate projective depths z.
  zEst = zeros(size(z));
  for j=1:nPts
    C = imPts(:,:,j);
    Ckron = zeros(3*nIm, nIm);
    for k = 1:nIm
      ei = zeros(nIm,1);
      ei(k) = 1;
      Ckron(:,k) = kron(ei, C(:,k));
    end

    % Minimize || Ckron * z - Mfac * Pfac(:,j) ||^2 subject to ||Ckron * z|| = 1
    Cz = Ckron * z(1, :, j)';
    b = Mfac * diag(Sfac) * Pfac(:,j);
    scl = norm(Cz).^-1;
    f0 = sum((scl * Cz - b).^2);
    
    if f0 > tol % Objective function is large enough, try to decrease it. 
      
      %%Require || Ckron z || = 1 = || diag(sCk) * vCk' * z_k ||
      % Set y = diag(sCk) * vCk' z_k, then ||y||_2 = 1.
      [uCk sCk vCk] = svd(Ckron, 0); sCk=diag(sCk);
      if DEBUG
        if sCk(1)/sCk(end) > 1000
          fprintf(2,' Cond number Ckron(%d): %e\n',j, sCk(1)/sCk(end));
        end
      end
      scl = sCk.^-1;
      
      % Minimize ||Pk * y - Mfac * diag(Sfac) * Pfac(:,j)|| with ||y||=1
      % where Pk is given by:
      Pk = Ckron * (vCk .* repmat(scl', nIm, 1));
      
      % Current guess is y0 = diag(sCk) * vCk' z_k.
      y0 = sCk.* (vCk' * z(1,:,j)');
      y0 = y0/norm(y0);
      
      % Simplify ||Pk * y - Mfac* diag(Sfac) * Pfac(:,j)|| to ||diag(sPk) * x - b ||,
      % subject to ||x|| = 1, where Pk = uPk * diag(sPk) * vPk', and vPk' * y = x
      [uPk sPk vPk] = svd(Pk, 0); sPk=diag(sPk);
      if DEBUG
        if sPk(1)/sPk(end) > 1000
          fprintf(2,' Cond number Pk(%d): %e\n',j, sPk(1)/sPk(end));
        end
      end
      b = uPk'*(Mfac * diag(Sfac) * Pfac(:,j));
      x0 = vPk' * y0;
      
      % Minimize ||diag(sPk) * x - b|| with ||x||=1, initial guess x0.
      % Compute the negative gradient direction
      delta = sPk .* (b - sPk .* x0);
      delta = delta - x0 * (x0' * delta);  % delta must be perp to x0.
      nDelta = norm(delta);
      delta = delta / max(nDelta, 1);
      
      % Try a step in the negative gradient direction of length alpha1
      alpha = 0;
      alpha1 = 1;
      cnt = 0;
      % Shrink alpha1 until objective function decreases.
      while alpha1 * nDelta > 10e-6
        x1 = x0 + alpha1 * delta;
        x1 = x1/norm(x1);
        f1 = sum((sPk .* x1 - b).^2);
        cnt = cnt+1;
        if f1 < f0 
          alpha = alpha1;
          break;
        end
        alpha1 = alpha1/2;
      end
      % Here alpha is either 0, or the objective function f1 is smaller
      % at this alpha and we can update the projective depths.
      
      % Compute projective depths for this alpha
      x = x0 + alpha * delta;
      x = x/norm(x);
      y = vPk * x;
      zEst(1,:,j) = vCk * (scl .* y);
      
    else  % Objective function small enough already.  Keep projective depths.
      zEst(1,:,j) = z(1,:,j);
    end
    
  end % End of z-estimation loop

  z = zEst;
  
end % End of projective factorization loop.
% Reshape projective depths from (1,nIm,nPts) to  (nIm, nPts).
z = reshape(z, nIm, nPts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Projective Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check out projective reconstruction.  
%%% Does it approximate the ground truth data up to a 3D Homography?
%%% The problem here is that we cannot expect to simply display projective
%%% reconstructions.  We first need to get the homography roughly right.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pest = Pfac;
Mest = Mfac * diag(Sfac(1:4));
MASSAGED_MEST = FALSE;

% Fit a homography from the estimated reconstruction, Pest, to the
% ground truth data P0.
[Hmassaged Sa] = linEstH3D(P0, Pest);

Pmassaged = Hmassaged * Pest;
XYZmassaged = Pmassaged./repmat(Pmassaged(4,:),4,1);

% Show 3D coords for transformed reconstruction.
% NOTE: The homography has be chosen using the ground truth.
% This just illustrates the consistency of the projective
% reconstruction with the original data.
figure(1); clf;
for k = 1:3
  subplot(1,3,k)
  hand1=plot(XYZmassaged(k,:), 'b', 'LineWidth', 2);
  set(gca,'FontSize', 14);
  hold on;
  hand2=plot(P0(k,:), 'r');
  if k == 2
    title('Recovered Projective Model (b), Ground Truth (r)');
  end
  xlabel('Point index');
  ylabel(sprintf('X(%d)\n',k));
end
fprintf(2,'Press any key to continue...');
pause;
fprintf(2,'ok\n');

% Summarize errors.
fprintf(2,'RMS Error: %f %f %f\n', sqrt(sum((P0(1:3,:) - XYZmassaged(1:3,:)).^2,2)/nPts));
mn0 = sum(P0(1:3,:),2)/nPts;
fprintf(2,'Data Point RMS variation: %f %f %f\n', sqrt(sum((P0(1:3,:) - repmat(mn0,1,nPts)).^2,2)/nPts));

%%%% Surface plot of transformed projective reconstruction.
figure(3); clf;
for k = 1:length(f)
  vf = XYZmassaged(:,f{k});
  patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
end
set(gca,'YDir', 'reverse');
axis vis3d; axis square; axis equal;
title('Projective Reconstruction (see NOTE)');
fprintf(2,'Rotate this figure.\n');
fprintf(2,'Press any key to continue...');
pause;
fprintf(2,'ok\n');

%%%% Surface plot of ground truth data
figure(10); clf;
for k = 1:length(f)
  vf = P0(:,f{k});
  patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
end
set(gca,'YDir', 'reverse');
axis vis3d; axis square; axis equal;
title('Ground Truth Data');
fprintf(2,'Rotate this figure.\n');
fprintf(2,'Press any key to continue...');
pause;
fprintf(2,'ok\n');

%%%% Show effect of various (non-affine) homographies
if HOMOG_DEMO
  for k = 1:10
    H = eye(4);
    H(4,1:3) = 0.02*(rand(1,3)-0.5);

    Ph = H * Pmassaged;

    %%%% Surface plot of data
    figure(4); clf;
    for k = 1:length(f)
      vf = Ph(1:3,f{k})./repmat(Ph(4,f{k}),3,1);
      patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
    end
    set(gca, 'FontSize', 14);
    axis vis3d; axis square; axis equal;
    title('Equivalent Projective Reconstruction');
    fprintf(2,'Press any key to continue...');
    pause;   fprintf(2,'ok\n');
  end
  clear H Ph;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finished checking out projective reconstruction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric Upgrade (and an attempt at self-calibration).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~MASSAGED_MEST  % Only do the following permutation of the Mest array once.
  if UNDO_RESCALE  
    % Undo the image transformation used for numerical stability, Knum
    % Not recommended.
    Mest = permute(reshape(Mest, [3, nIm, 4]), [1 3 2]);
    for kIm = 1:nIm
      Mest(:,:, kIm) = inv(reshape(Knum(:,:,kIm),3,3)) * reshape(Mest(:,:,kIm), 3, 4);
      Knum(:,:,kIm) = eye(3,3);
    end
    Mest = permute(Mest, [1 3 2]);
    Mest = reshape(Mest, [3*nIm, 4]);
  end
  Mest = reshape(Mest, [3, nIm, 4]);
  Mest = permute(Mest, [1 3 2]);
  MASSAGED_MEST = TRUE;
end

% Compute constraints on Qinf, the absolute dual quadric.
% Know Mest(i,:,k) * Qinf  * Mest(j,:,k)' = element i,j of K(:,:,k) * K(:,:,k)'
% K(:,:,k) the camera to image coordinate 3 x 3 matrix Knum * diag([f(k), f(k),1])
Acum = [];
bcum = [];
for iIm = 1:nIm
  
  % Compute constraints on Qinf from current image formation
  % matrix Mest(:,:,iIm).
  A = zeros(6,10);
  b = zeros(6,1);
  nc = 0;
  for i=1:3
    for j= i:3
      nc = nc+1;
      q = reshape(Mest(i,:,iIm),4,1) * reshape(Mest(j,:, iIm), 1,4);
      q = q + q' - diag(diag(q));
      A(nc,:) = [q(1,1:4) q(2,2:4) q(3,3:4) q(4,4)];
    end
  end
  KK = reshape(Knum(:,:,iIm),3,3) * diag([fMat(iIm) fMat(iIm) 1]);
  KK = KK * KK';
  % For (i,j)=(3,3) currently 6th row in A and b
  b = [KK(1, 1:3) KK(2, 2:3) KK(3,3)]';
  
  if SELF_CALIBRATE
    % For (i,j)=(1,1) and (2,2), 1st and 4th elements in A and b, get
    % combined constraint.  This cancels the dependence on the focal length.
    A(4,:) = A(4,:) - A(1,:);
    b(4) = b(4) - b(1);
    % Delete 1st row of both A and b to avoid extra constraint on (1,1) element.
    A = A(2:end, :);
    b = b(2:end);
  end
  
  % Concatenate constraints on Qinf across all images.
  if size(Acum,1)>0
    Acum = [Acum; A];
    bcum = [bcum; b];
  else
    Acum = A;
    bcum = b;
  end
  
end

% Solve A q = b for q and hence for Qinf.
[uA sA vA] = svd(Acum,0); sA = diag(sA);
% log10(sA(1:10)+eps)'
% Should be rank 10
if sA(10)/sA(1) < sqrt(eps)
  % Oh, oh.  Trouble!!!!
  fprintf(2,'Warning: System for Qinf nearly singular.\n');
  fprintf(2,'Press any key to continue...');
  pause;
  %  Solve A q = b, ignore singular element
  q = vA * (diag([sA(1:9).^-1 ; 0]) * (uA' * bcum));
  q0 = q;
  fprintf(2,'ok\n');
else
  %  Solve A q = b
  q = vA * (diag(sA.^-1) * (uA' * bcum));
  q0 = q;
end
if DEBUG
  fprintf(2, 'Residuals in solving for Qinf...');
  imStats(Acum * q - bcum);
end

% Search along the least constrained direction for Q for a rank 3 Q matrix.
CHATTY = TRUE;
foundMetric = FALSE;
foundZero = FALSE;
minEQ = nan;
if CHATTY
  fprintf(2,'Eigenvalues of estimated Qinf at each value of itNull:\n'); 
end
for itNull = -50:1:50
  q = q0 + itNull * vA(:,10);

  % Unpack q into symmetric 4x4 matrix
  Qinf = [q(1) q(2) q(3) q(4);
          q(2) q(5) q(6) q(7);
          q(3) q(6) q(8) q(9);
          q(4) q(7) q(9) q(10)];

  % Check out eigenvalues of estimated Qinf
  [uQ sQ vQ] = svd(Qinf); sQ = diag(sQ);
  sgn = diag(uQ' * vQ)';
  if CHATTY
    fprintf(2,'%f: ', itNull);
    fprintf(2,' %e',sQ(:).*sgn(:)); 
    fprintf(2,'\n');
  end
    
  if all(sgn(1:3)>0)
    if isnan(minEQ) 
      minEQ = sgn(4);
    elseif  minEQ * sgn(4) < 0 && ~foundZero
      if CHATTY
        fprintf(2, 'Detected zero crossing, itNull = %f\n',itNull);
      end
      foundMetric = TRUE;
      minEQ = sgn(3);
      itNull0 = itNull;
      foundZero = TRUE;
    end
  end
end

if ~foundMetric
  fprintf(2, 'Failed to do self-calibration\n');
else
  % Use itNull0
  q = q0 + itNull0 * vA(:,10);

  % Unpack q into symmetric 4x4 matrix
  Qinf = [q(1) q(2) q(3) q(4);
          q(2) q(5) q(6) q(7);
          q(3) q(6) q(8) q(9);
          q(4) q(7) q(9) q(10)];

  % Check out eigenvalues of estimated Qinf
  [uQ sQ vQ] = svd(Qinf); sQ = diag(sQ);
  sgn = diag(uQ' * vQ)';
  % Find closest non-negative definite rank 3 matrix to estimated Qinf
  kSv = 4; % Set smallest sv to set to zero
  sQinf = sQ;
  sQinf(kSv) = 0;
  Qinf = uQ * diag(sQinf) * uQ';  % Reconstructed rank 3 Qinf.

  %% Estimate homography for metric upgrade from this Qinf.
  sH = zeros(4,1);
  uH = zeros(4,4);
  n = 1;
  for k=1:4
    if k ~= kSv
      sH(n) = sQ(k); 
      uH(:,n) = uQ(:,k);
      n= n+1;
    end
  end
  uH(:,4) = uQ(:,kSv);
  sH(4) = 1;
  sH = sH .^ -0.5;
  H = (uH * diag(sH))';
  Hinv = uH' * diag(sH.^-1);

  if DEBUG
    fprintf(2, 'Sanity check; Transformed Qinf should be apprx diag([1 1 1 0])\n');
    H * Qinf * H'
  end

  %%% Given this estimated homography H, apply it the the projective
  %reconstruction.
  Hm2 = H;
  Pm2 = H * Pest;
  Pm2 = Pm2 ./ repmat(Pm2(4,:),4,1);

  %%%% Surface plot of metric reconstruction.
  figure(3); clf;
  for k = 1:length(f)
    vf = Pm2(:,f{k});
    patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
  end
  set(gca,'YDir', 'reverse', 'FontSize', 14);
  axis vis3d; axis square; axis equal;
  title('Euclidean Reconstruction');
  fprintf(2,' Rotate this reconstruction.\n');
  fprintf(2,...
          ' Note: The colour map is NOT supposed to be the same, just the shape.\n');

  pause;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Check out metric reconstruction.  
  %%% Does it approximate the ground truth data up to a 3D Homography?
  %%% The problem here is that we cannot expect to simply display projective
  %%% reconstructions.  We first need to get the homography roughly right.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Fit a rigid transform from the estimated reconstruction, Pest, to the
  % ground truth data P0.

  [Haff] = linEstAff3D(P0, Pm2);

  % FInd a nearby rigid transform;
  R0 = Haff(1:3,1:3);
  [Ur Sr Vr] = svd(R0); Sr = diag(Sr);
  Sr = sum(Sr)/3 * eye(3);
  R = Ur * Sr * Vr';
  Haff(1:3,:) = [R, R * inv(R0)* Haff(1:3,4)];

  Pmassaged = Haff* Pm2;
  XYZmassaged = Pmassaged./repmat(Pmassaged(4,:),4,1);

  % Show 3D coords for transformed reconstruction.
  % NOTE: The homography has be chosen using the ground truth.
  % This just illustrates the consistency of the projective
  % reconstruction with the original data.
  figure(1); clf;
  for k = 1:3
    subplot(1,3,k)
    hand1=plot(XYZmassaged(k,:), 'b', 'LineWidth', 2);
    set(gca,'FontSize', 14);
    hold on;
    hand2=plot(P0(k,:), 'r');
    if k == 2
      title('Recovered Euclidean Model (b), Ground Truth (r)');
    end
    xlabel('Point index');
    ylabel(sprintf('X(%d)\n',k));
  end
  fprintf(2,'Press any key to continue...');
  pause;
  fprintf(2,'ok\n');

  % Summarize errors.
  fprintf(2,'RMS Error: %f %f %f\n', sqrt(sum((P0(1:3,:) - XYZmassaged(1:3,:)).^2,2)/nPts));
  mn0 = sum(P0(1:3,:),2)/nPts;
  fprintf(2,'Data Point RMS variation: %f %f %f\n', sqrt(sum((P0(1:3,:) - repmat(mn0,1,nPts)).^2,2)/nPts));

  %%% Surface plot of metric reconstruction.
  figure(3); clf;
  for k = 1:length(f)
    vf = XYZmassaged(:,f{k});
    patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
  end
  set(gca,'YDir', 'reverse', 'FontSize', 14);
  axis vis3d; axis square; axis equal;
  title('Euclidean Reconstruction');
  fprintf(2,' Rotate this reconstruction.\n');
end


% Summary: For nIm perspective images, nIm >=3, we have demonstrated:
%   - Euclidean scene reconstruction from 3 or more orthographic images.
%   - The reconstruction of the viewing directions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End: projectiveMassageDino.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

