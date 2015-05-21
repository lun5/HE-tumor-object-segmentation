function [lines] = robustLineEstimator(p, theta, sigma, im, directedSeeds);
% [lines] = robustLineEstimator(p, theta, sigma, im, directedSeeds);
%
% Robust line fitting demo
% INPUT:
%  p Edgels are in variable p, an K by 2 matrix, with each row giving the
%    (x,y) location of the edgel.
%  theta  orientation of edgels in radians.  The edgel normal is
%         assumed to be in the direction (x,y) = (cos(theta), sin(theta)).
%         The normal is assumed to be directed along the image gradient
%         i.e. pointing towards the brighter side of the edge.
%  sigma  The spatial sigma used in the Canny edgel finder.
%  im     The original image (or cropped image) to be used for overlays.
%  directedSeeds (optional) A logical array indicating which edgels are
%         in oriented regions of the image.  Only these edgels will be
%         used as random seeds.  (Default: all edgels used, ie. directedSeeds
%         is all TRUE).
% OUTPUT:
%  lines  a 4xL array with the k-th column providing the endpoints,
%         namely,  (x0,y0, x1,y1)', for the k-th line found.

if ~exist('directedSeeds','var')
 directedSeeds = theta>=0.0;
end

lines = [];

FALSE = (0 == 1);
TRUE = ~FALSE;

DISPLAY_ON = FALSE;
CHATTY = FALSE;
echoPoints = 0;  % < 0 lineSupport only, >0 all Points + line support, = 0 none
projectWeights = TRUE;
useEdgeSign = FALSE;

Seeds = [1:size(p,1)]';
Nrml = [cos(theta) sin(theta)];
Px = p(:, 1);
Py = p(:, 2);
P = cat(2, Px, Py);

sigmaEst=2*sigma; 
sigmaX = 0.5 * sigma; alpha = 0.1;
rejectRadiusX = 4.0*sigmaX;
padX = 4.0;

nIts = 10; 
fitThres = 0.5;  fitRadius0 = 3.0; minWght0 = 1.25 * fitRadius0;
minOrientedNess0 = (1.0 - 0.1)/(1.0 + 0.1);
maxLine = 4000;
maxTries = 4 * maxLine;
sigmaDir = 0.15;  
plotIteration = nIts;
sigmaWin = 1.0;

sigmaWeights = 2;
streakWeights = 5;  % Odd
gFiltRad = round(2.5*sigmaWeights);
szWeightFilt = round(streakWeights + 2*gFiltRad+1);
if mod(szWeightFilt,2) == 0
  szWeightFilt = szWeightFilt +1;
end
tmp = zeros(1, szWeightFilt);
tmp((floor(szWeightFilt/2)-floor(streakWeights/2)):...
    (floor(szWeightFilt/2)+floor(streakWeights/2))) = 1.0;
gFilt = mkGaussian1(gFiltRad*2 + 1, sigmaWeights^2, gFiltRad+1, 'norm');
weightFilt = upConv(tmp,gFilt, 'zero');
weightFilt = weightFilt'/sum(weightFilt(:));

if DISPLAY_ON 
  hFig = figure(1);
  clf;
  image(im);
  colormap(gray(256));
  resizeImageFig(hFig, size(im), 2); hold on;
  hold on;
  if echoPoints > 0
    plot(Px,Py,'bo');
  end
  set(get(hFig,'CurrentAxes'),'Ydir','reverse');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over edgel selections, attempt to fit a line to each selection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kLine = 1; 
kTries = 1;
kUnSuccessful = 0; 
while(kLine <=maxLine & kTries <= maxTries & kUnSuccessful < 50)


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Generate an initial guess by sampling the edgels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Grab random edgel from directed seed edgels
  seedIndices = Seeds(directedSeeds);
  current = ceil(rand(1,1) * length(seedIndices));
  current = seedIndices(current);
 
  %% Form short segment covering this edgel, with endpoints p1, p2.
  p1 = P(current,:)';
  p2 = p1 + [Nrml(current, 2); -Nrml(current,1)];
  centerPt = (p1 + p2)/2.0;

  %% Determine tangent, and normal to this selected segment
  n = Nrml(current,:)';
  t = [n(2); -n(1)];
  c =  centerPt' * n; % distance from origin

  %% Debug output: Show current guess
  if DISPLAY_ON & plotIteration == 0
    projpoints=P*t;
    s1=min(projpoints); s2=max(projpoints);
    p1=c*n+s1*t; p2=c*n+s2*t;
    plot([p1(1),p2(1)],[p1(2),p2(2)],'y-');
    text(p1(1),p1(2),num2str(0));
    pause;
  end;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Select the edgels near the extended line n*x - c = 0
  %% This set will be updated, if necessary, as the line is
  %% re-estimated.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Get subset of edgels close to this current line.
  err = P * n - c;
  nPrev = n;
  cPrev = c;
  sumW = 0.0;
  nearLine = (abs(err) < rejectRadiusX + padX);
 
  %% Get the subset of these nearLine edgels close to the current segment.
  currentFitRadius = fitRadius0;  % Limit length of hypothesized segment.
  tangMid = centerPt' * t;
  tangSep = zeros(size(nearLine));
  tangSep(nearLine) = P(nearLine,:) * t - tangMid;
  nearSeg = nearLine & (abs(tangSep) <= currentFitRadius);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Robust fitting iterations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  foundLine = FALSE;
  for it=1:nIts
  
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% Select edgels near current guess.
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (it > 1)
      t = [n(2); -n(1)];
      tangMid = centerPt' * t;
      centerPt = c * n + tangMid * t;
       
      % Check if nearLine set needs to be recomputed.
      corners = [ -n -t; n -t; n t; -n t] * [rejectRadiusX; currentFitRadius];
      corners = reshape(corners, 2, 4) + centerPt * ones(1,4);
      recomputeNearLine = any(abs(nPrev' * corners - cPrev) >= rejectRadiusX+padX);

      % Recompute the near line set, if necessary
      if recomputeNearLine
        if CHATTY
          fprintf(2,'Recomputing near-line set...\n');
        end
        err = P * n - c;
        nPrev = n;
        cPrev = c;
        nearLine = (abs(err) < rejectRadiusX + padX);
      end
   
      % Compute near segment set.
      tangSep = zeros(size(nearLine));
      tangSep(nearLine) = P(nearLine,:) * t - tangMid;
      nearSeg = nearLine & (abs(tangSep) <= currentFitRadius);
    end
    %% Here nearSeg is true iff edgel is near the hypothesized segment

    %% Check for at least one close by segment.
    if ~any(nearSeg)
      if CHATTY
        fprintf(2, 'Hypothesized line received no support.\n');
      end
      foundLine = FALSE;
      break;
    end
   
    %% Set activeEdgelIndices to be the indices of the nearby edgels.
    tmp = 1:length(nearSeg);
    activeEdgelIndices = tmp(nearSeg);
    aP = P(activeEdgelIndices,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute errors and robust weights of active edgels 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute error residuals
    err = aP*n-c; 
    s = aP*t - tangMid;
    
    % Compute orientation residuals
    if useEdgeSign | (currentFitRadius <= 8.0)
      sameSign = (Nrml(activeEdgelIndices,:) * n) > 0;
    else
      sameSign = ones(length(activeEdgelIndices),1);
    end
    errOrient = sameSign.*(Nrml(activeEdgelIndices,:) * t)+(~sameSign);

    % Compute robust weights
    L = (err.^2)/(2*sigmaX^2) + (errOrient.^2)/(2*sigmaDir^2);
    maxL = ((rejectRadiusX/sigmaX)^2.0)/2.0;
    L(L<=maxL) = exp(-L(L<=maxL));
    L(L>maxL) = 0.0;
    W = L ./(alpha + L);

    % Enforce maximum radius of fitted segment
    W(abs(s)>currentFitRadius) = 0.0;

    % Sum of all weights contributing to current segment.
    prevSumW = sumW;
    sumW = sum(W);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check for sufficient weight and convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    converged = FALSE;
    if (it >= 2)
      if CHATTY
        fprintf(2, ' %d: sumW %f\n', kTries, sum(W));
      end
      if (sumW >= minWght0)
        foundLine = TRUE;
        if prevSumW > sumW - 0.5;
          % Converged
          converged = TRUE; 
        end
        prevSumW = sumW;
      else
        foundLine = FALSE;
        break;
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve weighted LSQ problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b=sum([W,W].* aP)';
    C=([W,W].* aP)'* aP /(sigmaX^2);
    sumW = max(eps, sum(W));
    C = C - b * (b'/(sigmaX^2 * sumW));
    Tang = [Nrml(activeEdgelIndices,2), -Nrml(activeEdgelIndices,1)];
    D =([W, W] .* Tang)' * Tang / (sigmaDir^2);
    [V E] = eig(C+D);
    E = diag(E);
    minIndex = (E == min(E));
    if (sum(minIndex)>1)
      minIndex(2) = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check for suitable seed (after 1 interation).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (it == 2)
      orientedNess = (E(~minIndex) - E(minIndex))/(eps + sum(E));
      if CHATTY
        fprintf(2, '    orientedness %e\n', orientedNess);
      end
      if (orientedNess >= minOrientedNess0)
        foundLine = TRUE; 
      else
        foundLine = FALSE; 
        break;
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute fitted line parameters, n, c and t.
    %% Switch sign of normal, if necessary, to match previous iteration.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = V(:, minIndex);
    if ( nPrev' * n < 0.0)
      n = -n;
    end
    t=[n(2);-n(1)];
    c = n' * b/sumW;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get extent of support for fitted line segment.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Project points and weights to line
    wIndices = W>0;
    projPoints = aP(wIndices, :)*t;
    projWeights = W(wIndices);

    % Find min and max of tangent projection, round to integers
    es1=min(projPoints); 
    es1 = floor(es1); 
    es2=max(projPoints);
    es2 = ceil(es2);
    dataExtent = es2-es1+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Form the 1D image, wIm, of edgel weights projected onto line segment.
    %% Use a Gaussian point spread function of standard dev sigmaWin.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigmaSqr = sigmaWin*sigmaWin;
    gFiltRad = round(3.0 * sigmaWin);
    x = -gFiltRad:1:gFiltRad;
    sclFilt = sum(exp(-(x .* x)/(2*sigmaSqr)));
    wIm = zeros(dataExtent + 2*gFiltRad, 1);

    % For each projected point with nonzero weight, add the edgel weight times
    % the point spread function to the weight image wIm.
    s0 = floor(projPoints);
    ds = projPoints - s0;
    for ks0 = 1:length(s0)
      % For the projected point with index ks0
      j0 = s0(ks0);
      w0 = projWeights(ks0);
      % Loop over the point spread function.
      for kf = (-gFiltRad):1:gFiltRad
        % Increment the weight image wIm.
        jIm = j0-es1+gFiltRad+1+kf; 
        wIm(jIm) = wIm(jIm) + w0 * exp( -(ds(ks0)+kf)^2/(2*sigmaSqr));
      end
    end
    % Rescale by the approximate integral of the point spread function.
    wIm = wIm/sclFilt;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Filter the 1D weight image, before looking for endpoints.
    %% Use the Gaussian * box filter, weightFilt computed above.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pad weight image with zeros, if necessary
    if length(wIm) < length(weightFilt)  
      % Pad wIm with zeros on end.
      wIm = [wIm(:) ; zeros(length(weightFilt) - length(wIm), 1)];
    end
    % Filter weight image
    fwIm = upConv(wIm, weightFilt, 'zero');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% A necessary condition for line support is that the filtered
    %% weight image must be larger than 0.3.  
    %% The left end of each of the line segments is then determined by the
    %% first point at which the weight image itself (i.e. wIm) becomes larger
    %% than 0.5 (so both fwIm>0.3 and wIm > 0.5 at the left endpoints). 
    %% Similarly the right end of each of the line segments is determined by
    %% the first point (from the right) at which the weight image itself 
    %% becomes larger than 0.5 (so both fwIm>0.3 and wIm > 0.5 at the
    %% right endpoints). 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    onIndicator = (fwIm > 0.3);
    lineOn = FALSE;
    % Erode left end of segment support intervals until wIm > 0.5.
    for k = 1:length(fwIm)
      if ~lineOn & onIndicator(k)
        lineOn = wIm(k) > 0.5;
        onIndicator(k) = lineOn;  %
      end
    end
    % Erode right end of segment support intervals until wIm > 0.5.
    lineOn = FALSE;
    for k = length(fwIm):-1:1
      if ~lineOn & onIndicator(k)
        lineOn = wIm(k) > 0.5;
        onIndicator(k) = lineOn;
      end
    end

    % Get segment beginning and end points only.
    startPts = onIndicator & [TRUE; ~onIndicator(1:(end-1))];
    stopPts =  onIndicator & [~onIndicator(2:end); TRUE];

    % Delete segments with only one wIm pixel of support.
    tmp = startPts & stopPts;
    onIndicator(tmp) = FALSE;
    startPts(tmp) = FALSE;
    stopPts(tmp) = FALSE;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Find contiguous segment of support containing (or closest to)
    %% the center point of the interval.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute integer array of difference between values of tangential
    % projection s, and the tangential projection of the current center point.
    sCurr = floor(centerPt' * t);
    sIm = (es1-1 - gFiltRad) - sCurr + [1:length(wIm)];
    % Get the index associated with the senter point.
    kCurr = sCurr - (es1-1 - gFiltRad);
    kCurr = max([kCurr, 1]);
    kCurr = min([kCurr, length(wIm)]);
    % If the center point is not within a supported line segment
    if ~onIndicator(kCurr)
      % Find nearest supported point, if any
      distLeft = length(wIm);
      for k = kCurr:-1:1
        if onIndicator(k)
          distLeft = kCurr - k;
          break;
        end
      end
      distRight = length(wIm);  
      for k = kCurr:length(wIm)
        if onIndicator(k)
          distRight = k-kCurr;
          break;
        end
      end
      % Move kCurr to nearest supported point
      if distLeft < distRight
        kCurr = kCurr - distLeft;
      elseif distRight < length(wIm)
        kCurr = kCurr + distRight;
      else
        kCurr = 0;  % No supported point found.
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Delete unconnected supported segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if kCurr == 0   % No supported pixels in wIm found.
      onIndicator(:) = FALSE;
    else
      onLine = TRUE;
      % Delete any extra subsegments on the right.
      for k = kCurr:length(wIm)
        if onLine & ~onIndicator(k)
          onLine = FALSE;
        elseif ~onLine & onIndicator(k)
          onIndicator(k) = FALSE;
        end
      end
      % Delete any extra subsegments on the left.
      onLine = TRUE;
      for k = kCurr:-1:1
        if onLine & ~onIndicator(k)
          onLine = FALSE;
        elseif ~onLine & onIndicator(k)
          onIndicator(k) = FALSE;
        end
      end
    end
    %% onIndicator contains exactly one supported segment, or none.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get beginning and end points of projected segment.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    startPts = onIndicator & [TRUE; ~onIndicator(1:(end-1))];
    stopPts =  onIndicator & [~onIndicator(2:end); TRUE];
    % Delete segments of length 1.
    tmp = startPts & stopPts;
    onIndicator(tmp) = FALSE;
    startPts(tmp) = FALSE;
    stopPts(tmp) = FALSE;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set minimum and maximum tangent projections
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(onIndicator)
      sMin = sIm(startPts) + sCurr;
      sMax = sIm(stopPts) + 1 + sCurr;
    else
      sMin = nan;
      sMax = nan;
      foundLine = FALSE;
      break;
    end

 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot projected weight image and endpoint determination results:
    %%   projected weight image wIm (k)
    %%   filtered wIm (b)
    %%   selected support segment (y)
    %%   start and stop points (g,r)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DISPLAY_ON & ((it >= plotIteration) | (converged & plotIteration <= nIts))
      figure(2); clf; 
      plot(sIm, wIm, 'k');
      hold on;
      plot(sIm, fwIm, 'b');
      plot(sIm, onIndicator, 'y');
      plot(sIm, startPts, 'g');
      plot(sIm, stopPts, 'r');
      axis([min(sIm), max(sIm), 0, 1.1]);
      title('Projected edgel weights and selected endpoints(r,g)');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Show fitted extended line and/or fitted line segment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  DISPLAY_ON & (converged | (it >= plotIteration))
      pt = centerPt' -n' * (centerPt' * n  - c);
      sCurr = centerPt' * t;
      mx = [size(im,2); size(im,1)];
      endPts = cropLineInBox(n, -c, [0 0 mx']);
      if ~isnan(endPts(1,1))
        figure(hFig); hold on;
        % Covering line
        plot(endPts(:,1), endPts(:,2),'y-', 'LineWidth', 0.5);
        % Plot segment
        qw = [pt + (sMax - sCurr) * t'; pt+(sMin-sCurr) *t'];
        plot(qw(:,1), qw(:,2),'g-', 'LineWidth', 2.0);
        if plotIteration < nIts
          text(qw(1,1),qw(1,2),num2str(it));
        end
        drawnow;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Prepare for next iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tangMid = (sMin + sMax)/2.0;
    centerPt = c * n + tangMid * t;
    currentFitRadius = (sMax - sMin)/2.0;
    
    if converged
      break;
    end

    %% Extend line segment support for next iteration.
    if it < nIts
      currentFitRadius = max([currentFitRadius + fitRadius0 * 2,...
                              currentFitRadius * 2]);
    end
  
  end  % End of robust estimation iteration loop

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Finish processing current random seed.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (foundLine)

    if ~isnan(sMin) & ~isnan(sMax)
      s = aP * t;
      W(s<sMin | s>sMax) = 0.0;
      tangMid = (sMin + sMax)/2.0;
      centerPt = c * n + tangMid * t;
      currentFitRadius = (sMax - sMin)/2.0;
      seg =repmat(centerPt, 1,2) + [t t]*diag([sMax-tangMid; sMin-tangMid]);
      lines = [lines seg(:)];
    end

    if DISPLAY_ON & echoPoints
      hFig = figure(1);
      hold on;
      plot(aP(W>=fitThres,1),aP(W>=fitThres,2),'ro');
      drawnow;
    end


    %% Remove covered edgels from seeds
    fitEdgels = zeros(size(Seeds));
    fitEdgels(activeEdgelIndices) = (W>=fitThres); 
  
    Seeds = Seeds(~fitEdgels);
    P = P(Seeds, :); Px = P(:, 1); Py = P(:, 2);
    Nrml = Nrml(Seeds,:);
    directedSeeds = directedSeeds(Seeds);
    Seeds = 1:size(P,1);
    kUnSuccessful = 0;
    kLine = kLine+1;

  else

    % Keep track of number of successive unsuccessful guesses 
    kUnSuccessful = kUnSuccessful + 1;

  end

  % Total number of seeds tried;
  kTries = kTries+1;

end  % over guesses.

kLine = size(lines,2);
display(kLine)
