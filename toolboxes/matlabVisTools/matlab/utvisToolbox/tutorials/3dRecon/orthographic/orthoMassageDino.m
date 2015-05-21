%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script file: orthoMassageDino.m 
%
% Demonstrate 3D affine and Euclidean reconstructions from corresponding points in
% orthographic projection. 
% 
% ADJ, Nov. 03. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear
global TRUE FALSE;
global matlabVisRoot;

TRUE = 1==1;
FALSE = ~TRUE;

if isempty(matlabVisRoot)
  dir = pwd;
  cd '/h/51/jepson/pub/matlab'; % CHANGE THIS
  %cd 'C:/work/Matlab'
  startup;
  cd(dir);
end
reconRoot = [matlabVisRoot '/utvisToolbox/tutorials'];

addpath([reconRoot '/3dRecon/utils']);
addpath([reconRoot '/3dRecon/orthographic']);

% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

% Parameters
DEBUG = FALSE;
sigma = 2.0;         % Std Dev of noise (in pixels) in point locations
nIm = 3;             % Number of data images to use (try 2-10).
                     % nIm = 2 is enough for Affine reconstruction, but not
                     % for Euclidean reconstruction.
NeckerReversal = FALSE;   % Flip reconstruction in depth (the sign is ambiguous).

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

%%% Surface plot of data.
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
dMat = [50  -50    0   50   -40  -50   10   20     2.5    10
         0    0 -100  -20    10    0  -50  -20    10     -10
      -150 -150 -145 -160  -145  155 -150 -160  -140    -145];
% Focal lengths for each of these 10 cameras.
fMat = [1.25 1.25 0.75 1.25 1.25 1.20 0.8 1.0 1.20 1.0];
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
  K = diag([f1, f1, 0]);  % Scaled-Orthographic
  M = K * [R -R*d];
  if size(M_GT,1)>0
    M_GT = [M_GT; M(1:2, 1:3)];
  else
    M_GT = M(1:2, 1:3);
  end

  % Compute the orthographic image locations
  p = M * P0;
  
  % Add imaging noise.
  p(1:2,:) = p(1:2,:) + sigma * randn(2,nPts);

  % Concatenate image points over multiple frames.
  if size(imPts,1)>0
    imPts = cat(3, imPts, p(1:2,:));
  else
    imPts = p(1:2,:);
  end

  % Show current noisy image
  figure(10); clf;
  imStats(p(1,:));
  imStats(p(2,:));
  h = plot(p(1,:), p(2,:), '.b');
  set(gca,'YDir', 'reverse', 'FontSize', 14);
  axis([-200 200 -200 200]);
  title(sprintf('Image %d',k));
  fprintf(2,'Press any key to continue...');
  pause; fprintf(2,'ok\n');
  
end  % Forming images.

% Reorder the matrices imPts in the form nIm by nPts
if size(imPts,2) == nPts && size(imPts,3) == nIm
  imPts = permute(imPts,[1, 3, 2]);
end
imPtsSave = imPts;  % Save the image points in case we mess up...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END of data generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do Orthographic Factorization: Recover Affine Shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imPts = imPtsSave;

% Do factorization 
D = reshape(imPts, 2 * nIm, nPts);
if size(D,2)<=size(D,1)
  [U S V] = svd(D,0);  S = diag(S);
else
  [V S U] = svd(D',0);  S = diag(S);
end

% Print size of 4th singular value relative to first.  Is exactly zero without
% noise and round-off errors. 
fprintf(2, 'sv(4)/sv(1) %e\n', S(4)/S(1));

% Plot singular values.  Data matrix should be essentially rank 3, if
% imaging noise is small enough.
figure(1); clf; plot(S, '-*', 'MarkerSize', 14, 'LineWidth', 2);
set(gca, 'FontSize', 14);
title('Singular Values of Data Matrix');
xlabel('Singular Value Index');
ylabel('Singlur Value');
pause(0.1);

% Extract estimates for affine shape.
Mest = U(:,1:3) * diag(S(1:3));
Pest = V(:,1:3)';

% Show estimate for affine shape.
figure(1); clf; 
showWire(Pest', f);
set(gca, 'FontSize', 14);
xlabel('X'); ylabel('Y'), zlabel('Z');
title('Affine Estimation of Dino');
pause(0.1);

%%% Surface plot of data.
figure(3); clf;
for k = 1:length(f)
  vf = Pest(:, f{k});
  patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
end
set(gca,'YDir', 'reverse', 'FontSize', 14);
axis vis3d; axis square; axis equal;
title('Affine Estimation of Dino');
xlabel('X'); ylabel('Y'), zlabel('Z');
fprintf(2,'Rotate this figure.\n');
fprintf(2,'Press any key to continue...');
pause; fp rintf(2,'ok\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try Euclidean Reconstruction.
% For nIm == 2 we will get the bas relief ambiguity (rotation versus depth).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Square pixels:
%  Seek a matrix Q such that the rows of MestQ = Mest *Q for image k have the same
%  length.  That is norm(MestQ(2*k-1,:)) = norm(MestQ(2*k, :)) for k = 1..nIm.
%  We also require these rows to be orthogonal, i.e. 
%  MestQ(2*k-1,:) * MestQ(2*k, :) = 0, for k = 1..nIm.
%  These constraints ensure the pixels are square in each image.
% Overall scale factor:
%  Finally to resolve the overall scale (which is ambiguous), we will
%  set the norm of the rows of MestQ for the first image to 1.  That is,
%  norm(MestQ(1,:)) = norm(MestQ(2,:)) = 1.

% These constraints are encoded in terms of the values of particular
% rows  and columns in the matrix:
%   (Mest * Q) (Mest * Q)') = Mest * Q * Q' * Mest'
% We Set QQT = Q * Q'.  It turns out we need to solve a linear system
% for the elements of QQT.  Note QQT is a symmetric 3x3 matrix, so we
% only need to recover 6 coefficients.  These will be stored in the
% q2, where:
%   QQT = [q2(1) q2(2) q2(3); q2(2) q2(4) q2(5); q2(3) q2(5) q2(6)]
% We build the linear system for q2 next, in the form:  C * q2 = b;

C = zeros(2*(nIm-1)+3, 6);
b = zeros(2*(nIm-1)+3,1);
kC = 0;  % The constraint counter.
for kIm = 1:nIm
  kM = 2*(kIm-1);
  % Get rows of Mest for image kIm.
  Mk = Mest(kM+(1:2),:);
  % Build the constraints on the length of the rows of Mk.
  if kIm == 1
    % For the first image, require the length of both rows is 1.
    % This sets the overall scale factor.
    for r = 1:2
      outer = Mk(r,:)' * Mk(r,:);
      outer = outer + outer' - diag(diag(outer));
      C(kC+r,:) = [outer(1,1:3) outer(2,2:3) outer(3,3)];
      b(kC+r) = 1;
    end
    kC = kC+2;
  else  % For images 2 through nIm, require the lengths of both rows of
        % Mk are equal.  This allows the scale of each image to be different.
    for r = 1:2
      outer = Mk(r,:)' * Mk(r,:);
      outer = outer + outer' - diag(diag(outer));
      if r==1
        C(kC+1,:) = [outer(1,1:3) outer(2,2:3) outer(3,3)];
      else
        C(kC+1,:) = C(kC+1,:) - [outer(1,1:3) outer(2,2:3) outer(3,3)];
      end
    end
    b(kC+1) = 0;
    kC = kC+1;
  end
  % Make sure the two rows of Mk are orthogonal.
  outer = Mk(1,:)' * Mk(2,:);
  outer = outer + outer' - diag(diag(outer));
  C(kC+1,:) = [outer(1,1:3) outer(2,2:3) outer(3,3)];
  b(kC+1) = 0;
  kC = kC+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Q.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nIm>2
  
  % If the camera locations are not (nearly) coplanar, we should be able to solve uniquely
  % for Q. 
  [Uc Sc Vc] = svd(C,0); Sc = diag(Sc);
  log10(Sc)  % Is C singular?  It will be if the cameras are coplanar
             % (and there is no noise).
  tmp = Uc' * b;
  tmp = tmp./Sc;
  q2 = Vc * tmp;
  log10(abs(C*q2 - b))  % Did we solve the linear system?
  
  % Build the symmetric matrix QQT def= Q * Q' from q2.
  QQT = [q2(1:3)'; q2(2) q2(4:5)'; q2(3) q2(5) q2(6)];
  
  % Factor QQT into Q and Q', by first finding the SVD.
  [Uq Sq Vq] = svd(QQT,0); Sq = diag(Sq)
  sgn = diag(Vq' * Uq)  % These MUST all be positive, otherwise we get a
                        % complex-valued factor.
  
  % Panic if the signs don't work out.
  if ~all(sgn>0)
    error('Panic, QQT not postive definite');
  else
    % Get the factors Q and Q' from the SVD of QQT.
    Q = Uq * diag(Sq.^0.5) * Vq';
    if NeckerReversal  % Depth reversal
      Q = -Q;          % Note the sign of Q is ambiguous.  Changing the sign
    end                % of Q effectively reverses the reconstruction in depth.
    MestQ = Mest * Q;
    sqrt(sum((MestQ).^2,2))  % This should be constant, in pairs.  The
                             % constants should be the scale factors in fMat.
    PestQ = inv(Q) * Pest;
    
    % Show the Euclidean reconstruction (shape is known up to a 3D
    % similarity transform, i.e. up to an unknown global rotation, translation,
    % and scale change.
    figure(1); clf; 
    showWire(PestQ', f);
    xlabel('X'); ylabel('Y'), zlabel('Z');
    title('Euclidean Reconstruction');
    
    fprintf(2,'Press any key to continue...');
    pause; fprintf(2,'ok\n');
  end
  
elseif nIm==2
  
  % We cannot solve for Q, we have only 5 constraints for 6 unknowns in QQT.
  % This is the bas relief ambiguity in which the overall depth
  % variation of the object is linked to an unkown rotation.
  
  % Since C is only 5 x 6 we only have 5 singular values.
  [Uc Sc Vc] = svd(C); Sc = diag(Sc);
  log10(Sc)  % The 5 singular values are nonzero...unless the two images are
             % have the same projection direction.
  
  % Let's solve C*q2 = b in the least squares sense.
  tmp = Uc' * b;
  tmp = tmp(1:5)./Sc(1:5);
  q2part = Vc(:, 1:5) * tmp;  % Least squares solution.
  % Since C is 5 x 6, we have C * (q2part + alpha * q2null) = b is also
  % a least squares solution for any real number alpha, where q2null is
  % a right null vector for C, i.e. C * q2null = 0.
  q2null = Vc(:,6);  
  
  % How well have we solved the system for q2?
  log10(abs(C*q2part - b))

  % Find a range of values for alpha for which QQT has real factors Q
  % and Q'.
  
  % Try alpha=0.
  alpha = 0;
  q2 = q2part+alpha*q2null;
  QQT = [q2(1:3)'; q2(2) q2(4:5)'; q2(3) q2(5) q2(6)];
  [Uq Sq Vq] = svd(QQT,0); Sq = diag(Sq);
  sgn = diag(Vq' * Uq);   % We have real factors iff these are all non-negative.
  sgnOkZero = all(sgn>0);
  if sgnOkZero
    % We can factor the case when alpha=0.  Let's do it. 
    Q = Uq * diag(Sq.^0.5) * Vq';
    if NeckerReversal
      Q = -Q;
    end
    PestQ = inv(Q) * Pest;
    % Keep track of the aspect ratio of the reconstruction.
    [Vp Sp Up] = svd(PestQ', 0); Sp = diag(Sp);
    zAspect = Sp(3)/Sp(1);  % This is actually the inverse aspect ratio,
                            %to avoid divide by zeros.
  end

  % Check positive and negative values of alpha for ranges in
  % which q2 = q2part+alpha*q2null provides a matrix QQT with real factors.
  pAlpha = [];
  nAlpha = [];
  nAspect = [];
  pAspect = [];
  for la10 = -10:0.1:10
    alpha = 10.0^la10;  % Use equal spacing in log10(alpha) to check a large range of values
    
    % Check out factorization of corresponding QQT for this value of alpha.
    q2 = q2part+alpha*q2null;
    QQT = [q2(1:3)'; q2(2) q2(4:5)'; q2(3) q2(5) q2(6)];
    [Uq Sq Vq] = svd(QQT,0); Sq = diag(Sq);
    sgn = diag(Vq' * Uq);
    if all(sgn>0) % Factorization of QQT will be ok.
      % Remember the value of alpha.
      pAlpha = [pAlpha alpha];
      % Do the factorization (because we want to compute the aspect
      % of the reconstruction below).
      Q = Uq * diag(Sq.^0.5) * Vq';
      if NeckerReversal
        Q = -Q;
      end
      % Make note of the (inverse) aspect ratio of the reconstruction.
      PestQ = inv(Q) * Pest;
      [Vp Sp Up] = svd(PestQ', 0); Sp = diag(Sp);
      pAspect = [pAspect Sp(3)/Sp(1)];
    end
    
    % Do the same for alpha with the opposite sign (-alpha).
    alpha = -alpha;
    q2 = q2part+alpha*q2null;
    QQT = [q2(1:3)'; q2(2) q2(4:5)'; q2(3) q2(5) q2(6)];
    [Uq Sq Vq] = svd(QQT,0); Sq = diag(Sq);
    sgn = diag(Vq' * Uq);
    if all(sgn>0) % Factorization of QQT will be ok.
      % Remember the value of alpha.
      nAlpha = [alpha nAlpha];
      % Do the factorization (because we want to compute the aspect
      % of the reconstruction below).
      Q = Uq * diag(Sq.^0.5) * Vq';
      if NeckerReversal
        Q = -Q;
      end
      % Make note of the (inverse) aspect ratio of the reconstruction.
       PestQ = inv(Q) * Pest;
      [Vp Sp Up] = svd(PestQ', 0); Sp = diag(Sp);
      nAspect = [ Sp(3)/Sp(1) nAspect];
    end
  end  % Finished checking out a range of positive and negative alpha's.
  
  % Put together a range of alpha's for which QQT has real factors.
  if sgnOkZero 
    alpha = [0 pAlpha];
    aspect = [zAspect pAspect];
    if size(nAlpha,1)>0
      alpha = [nAlpha alpha];
      aspect = [nAspect aspect];
    end
  elseif size(nAlpha,1)==0
    alpha = pAlpha;
    aspect = pAspect;
  elseif size(pAlpha,1) == 0
    alpha = nAlpha;
    aspect = nAspect;
  else
    fprintf(2, 'Found two seqments of alpha\n');
    alpha = [nAlpha pAlpha];
    aspect = [nAspect pAspect];
  end
  
  % Let's look at some of the reconstructions for different alpha's
  if length(alpha) > 0
    
    % Select just a few values of alpha to display.
    k = ceil(length(alpha)/5);
    kal = 1:k:length(alpha);
    if (length(alpha)-kal(end))/k > 0.5
      kal = [kal length(alpha)];
    else
      kal(end) = length(alpha);
    end
    
    % Display the affine reconstructions for the selected alpha's.
    for alp = alpha(kal)
      % Do the factorization.
      q2 = q2part+alp*q2null;
      QQT = [q2(1:3)'; q2(2) q2(4:5)'; q2(3) q2(5) q2(6)];
      [Uq Sq Vq] = svd(QQT,0); Sq = diag(Sq);
      sgn = diag(Vq' * Uq);  % Should be all positive (we don't check).
      Q = Uq * diag(Sq.^0.5) * Vq';
      if NeckerReversal
        Q = -Q;
      end
      % Find the camera and shape matrices.
      MestQ = Mest * Q;
      PestQ = inv(Q) * Pest;
      
      % Show the affine reconstruction. 
      % Note the shape is elongated or flattened, depending on the value
      % of alpha.  This is the bas relief ambiguity, where the amount of
      % rotation in depth is confounded with the amount of depth variation.
      %%% Surface plot of data.
      figure(1); clf;
      for k = 1:length(f)
        vf = PestQ(:, f{k});
        patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
      end
      set(gca,'YDir', 'reverse', 'FontSize', 14);
      axis vis3d; axis square; axis equal;
      title(sprintf('Affine Reconstruction (alpha = %9.2e)', alp));
      xlabel('X'); ylabel('Y'), zlabel('Z');
      fprintf(2,'Rotate this figure.\n');
      
      if FALSE
      figure(1); clf; 
      showWire(PestQ', f);
      xlabel('X'); ylabel('Y'), zlabel('Z');
      title(sprintf('Affine Reconstruction (alpha = %9.2e)', alp));
      end
      
      fprintf(2,'Press any key to continue...');
      pause; fprintf(2,'ok\n');
    end
  end
end  % Done solving for Q.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show surface plot of reconstruction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nIm ==2
  % For nIm == 2, choose the minimum aspect ratio depth estimate.  This
  % avoids reconstructions that are flat, or stretched in 1D due to the
  % bas relief ambiguity we just saw.
  [asp k] = max(aspect);
  alp = alpha(k);

  % Factor QQT, solve for Q
  q2 = q2part+alp*q2null;
  QQT = [q2(1:3)'; q2(2) q2(4:5)'; q2(3) q2(5) q2(6)];
  [Uq Sq Vq] = svd(QQT,0); Sq = diag(Sq)
  sgn = diag(Vq' * Uq)
  Q = Uq * diag(Sq.^0.5) * Vq';
  if NeckerReversal
    Q = -Q;
  end
  % Generate camera and shape matrices.
  MestQ = Mest * Q;
  sum((MestQ).^2,2)
  PestQ = inv(Q) * Pest;
  
  % Show wire-frame reconstruction.
  figure(1); clf;
  showWire(PestQ', f);
  title('Selected Affine Reconstruction');
  xlabel('X'); ylabel('Y'), zlabel('Z');
  fprintf(2,'Press any key to continue...');
  pause; fprintf(2,'ok\n');
  
  %%% Surface plot of selected affine reconstruction.
  figure(3); clf;
  for k = 1:length(f)
    vf = PestQ(:,f{k});
    patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
  end
  set(gca,'YDir', 'reverse');
  axis vis3d; axis square; axis equal;
  xlabel('X'); ylabel('Y'), zlabel('Z');
  title('Selected Affine Reconstruction');
  fprintf(2,'Rotate this figure.\n');
  
  fprintf(2,'Press any key to continue...');
  pause; fprintf(2,'ok\n');
  
elseif nIm > 2  % Show Euclidean reconstruction.
  
  % We have already solved for Q, MestQ, and PestQ above. 
  % Q is unique up to the sign (assuming the projection directions are
  % not coplanar).  Changing the sign of Q flips the reconstruction in depth.
  % This depth reversal is called the Necker ambiguity.
  
  %%% Surface plot of Euclidean reconstruction.
  figure(3); clf;
  for k = 1:length(f)
    vf = PestQ(:,f{k});
    patch(vf(1,:)', vf(2, :)', vf(3,:)', -vf(3,:)');
  end
  set(gca,'YDir', 'reverse', 'FontSize', 14);
  axis vis3d; axis square; axis equal;
  xlabel('X'); ylabel('Y'), zlabel('Z');
  title('Euclidean Reconstruction');
  fprintf(2,'Rotate this figure.\n');
  
  fprintf(2,'Press any key to continue...');
  pause; fprintf(2,'ok\n');
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check out the focal length estimate (scale factor).
% Does it approximate the ground truth data?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMest = sqrt(sum(reshape(sum((MestQ).^2,2), 2, nIm), 1)/2);
% Match over-all scale factor for comparison purposes.
% We couldn't do this without some absolute scale factor.  In that 
% case we would still recover the relative scales of all the images.
sclFac = (fMest * fMat(1:nIm)')/(fMest * fMest')
figure(4); clf;
plot(sclFac*fMest, 'b*-'); hold on;
plot(fMat(1:nIm), 'ro-');
plot(100*abs(sclFac*fMest-fMat(1:nIm))./fMat(1:nIm), 'k');
set(gca,'FontSize', 14);
title('Scale estimation (Est(b) True(r) %Error(k))');
xlabel('Image Number');
ylabel('Scale or % Error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check out affine reconstruction.  
% Does it approximate the ground truth data up to a 3D affine transform?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit a 3d affine transform for the estimated reconstruction, PestQ, from the
% ground truth data P0.  That is, solve (A b) * P0  = PestQ

Ab = ( inv(P0 * P0') * (P0 * PestQ'))';
Paff = Ab*P0;

% Show the error (in units of pixel size in X,Y and ALSO Z)
errPnts = sqrt(sum((Paff - PestQ).^2 ,1));
figure(1); clf; plot(errPnts);
title(sprintf('Affine Error in Recovered Points'));
xlabel('Point index');
ylabel('Euclidean Distance (pixels)');
pause(0.1);

% Show individual coordinates of the recovered shape, X,Y and Z.
% Compare these with the ground truth.
figure(2); clf;
coord = 'XYZ';
for k = 1:3
  subplot(1,3,k);
  plot(PestQ(k,:),'b', 'MarkerSize', 14, 'LineWidth', 2);
  hold on;
  plot(Paff(k,:),'r', 'MarkerSize', 14, 'LineWidth', 2);
  set(gca, 'FontSize', 12);
  title(sprintf('Recovered Coord (b), True(r)'));
  xlabel('Point index');
  ylabel(['3D coord ' coord(k)]);
end
pause(0.1);

%% What about the Euclidean reconstruction?
% The above error estimates compare the affine reconstruction to
% the ground truth. To check the Euclidean reconstruction we should find
% a similarity transform of P0  (i.e. (sR d) * P0, s>0 a scale
% factor, R a rotation matrix, and d a translation vector) for which
% (sR d) P0 provides the best approximation of PestQ.
% This is a nonlinear optimization problem (because of the rotation
% R).  We do not do this here.
%
% Instead, we simply check how close to a rigid transform is the best affine
% transform from P0 to PestQ.

% Find the svd of the A portion of the affine transform given
% by Ab = (A b) (a 3 x 4 matrix)
[Ua Sa Va] = svd(Ab(:, 1:3)); Sa = diag(Sa);
Sa  % If we have a similarity transform, then the singular values will
    % all be equal.
Sa(3)/Sa(1)  % This should be near 1 if we have a decent Euclidean
             % reconstruction... i.e A approx= sR for some rotation
             % matrix R.  It will probably not be near 1 when nIm == 2
             % because of the bas relief ambiguity.
             
% A guess at the rotation/reflection matrix is obtained by setting all
% the singular values of A to be equal, providing
R_est = Ua * Va';
             
% Finally, we need the determinant of R_est to be positive. Otherwise we
% have reconstructed the depth reversed scene (which is just the Necker
% ambiguity again).
det(R_est)  % Will be +-1 since the determinant of the product of two
            % matrices is the product of the determinants 
            % (i.e. det(M N) = det(M)det(N)) and the determinant of an
            % orthogonal matrix is +-1 (i.e. det(Ua) = det(Va) = +-1 ). 

if det(R_est) < 0
  NeckerReversal = ~NeckerReversal;
  % In order to reverse the depth, paste in from comment 
  % "Solve for Q" above, all the way to the end here.
  % (Don't simply rerun the whole file, since the code will pick new
  % noise values, and this may also cause the estimate to reverse in depth.)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check out the recovered projection directions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errProjDir = zeros(1, nIm);
for kIm = 1:nIm
  
  % Estimate the projection direction for the k-th image.
  % Unpack the k-th camera matrix
  Mk = MestQ((2*kIm-1):(2*kIm),:);
  % Find the right null vector of Mk.
  [Vm Sm Um] = svd(Mk'); Sm = diag(Sm);
  % Vm(:,3) is the recovered projection direction of the k-th camera
  projDirEst = Vm(:,3);
  
  % Do the same for the ground truth imaging matrices
  % Unpack the k-th camera matrix
  Mk = M_GT((2*kIm-1):(2*kIm),:);
  % Find the right null vector of Mk.
  [Vm Sm Um] = svd(Mk'); Sm = diag(Sm);
  % Vm(:,3) is the recovered projection direction of the k-th camera
  projDir_GT = Vm(:,3);
  
  % The estimated and true directions should be related by the rotation
  % in the Euclidean transformation.
  errProjDir(kIm) = asin(norm(cross(projDirEst, R_est * projDir_GT)))*180/pi;
end

figure(5); clf;
plot(errProjDir, '*-b', 'MarkerSize', 14, 'LineWidth', 2);
set(gca,'FontSize', 14);
title('Error in estimated projection directions');
xlabel('Image Number');
ylabel('Error in Degress');
ax = axis;
axis([ax(1:2) 0 ax(4)]);

% Summary: For nIm orthographic images, nIm >=2, we have demonstrated:
%   - Euclidean scene reconstruction from 3 or more orthographic images.
%   - Affine scene reconstruction from nIm = 2 orthographic images, and the
%     associated bas relief ambiguity.
%   - The Necker (depth reversal) ambiguity.
%   - The reconstruction of the viewing directions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End: orthoMassageDino.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
