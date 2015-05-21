%%%   Setup for Motion Tutorial
%%%  File: motionTutorial.m
%%%  Matlab script file
%%%  Date: Oct, 03
%%%
%%% Compute brightness constancy constraints, and estimate
%%% translational motion.  Demonstrate results by tracking a
%%% small image region

%%% Dependencies
%%%   Image sequences:
%%%      matlabVisRoot/
%%%         images/seq/fleetface1/ 
%%%         images/seq/dudekface1/
%%%   Toolboxes:
%%%      iseToolbox/
%%%      utvisToolbox/file/
%%%      utvisToolbox/tutorials/motionTut/

%%% Things still to do:
%%%   Similarity or affine warps.
%%%   Robust estimation of motion.
%%%   Coarse to fine motion estimation.

%%% Author: ADJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Check Path and Constants  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

global matlabVisRoot;
if ~exist('matlabVisRoot', 'var') || isempty(matlabVisRoot)
  %% Set matlabVisRoot to the directory containing utvisToolbox and iseToolbox.
  matlabVisRoot='/h/51/jepson/pub/matlab/';
end

%  Check path.  If these aren't on your path, then you will need
%  to add these toolboxes to your path.  See ~jepson/pub/matlab/startup.m
which resizeImageFig  %% should be utvisRoot/file/resizeImageFig.m
which showIm          %% should be iseToolbox/pyrTools/showIm.m

addpath([matlabVisRoot '/utvisToolbox/tutorials/motionTut/']);
which cropImage;  %% should be utvisToolbox/tutorials/motionTut/cropImage.m

TRUE = 1 == 1;
FALSE = ~TRUE;

%%% Default Parameters

SEQ_PREVIEW = FALSE;  %(default) FALSE 
% Set to TRUE if you want to preview
% the first few frames of the sequence 

chooseSeq = 1;  % Use fleetface (1) or dudekface sequence (2)

PLOT_BCC = TRUE; %(default) TRUE; 
% Plot the brightness constancy constraints.
% The plotting slows things down, you may wish to turn it off to speed
% things up.  But only do this after you have understood the plots.

nPlotBCC = 200; % (default) 200;
% Number of brightness constancy constraints to plot.  These
% constraints are randomly sampled.  Too many constraints makes
% the plots un-interpretable (i.e. all black).

ZERO_PREDICT = TRUE; %(default) TRUE;  
% For each new frame we need to predict the location of the matching
% patch.  When ZERO_PREDICT is TRUE we simply use the aame position as
% in the previous frame  When ZERO_PREDICT is FALSE we use the position 
% predicted from both the previous frame and the previous speed.  Start
% with this set to TRUE, as this provides worse predictions and therefore
% exercises the convergence properties of our motion estimation algorithm
% more fully.

USE_MOUSE = FALSE; % (default) FALSE
% Set to true if you wish to mouse in a box to
% determine the starting region in the first frame.  Start with this
% being FALSE, to see the examples we have set up with you.

mxRewarpIts = 3;  % (default) 3
% This specifies the max number of iterative rewarp steps.
% Try the default value 3.  Compare this with mxRewarpIts = 1 to
% see what happens with no iterative rewarping (just the linear
% brightness constancy constraints).  The effect is largest with
% smaller values of sigmaBlur (eg. 1 or 2) and ZERO_PREDICT = TRUE (
% so the initial guesses are worse.

sigmaBlur = 2;  % (default) 2
% Sigma for all gaussian and gaussian derivative filters.
% You can try larger (sigmaBlur = 4, 8) or smaller eg sigmaBlur=0.5, 1
% values of this.  Remember the linearized brighness constancy
% constraints approximte the image as linear over a distance of roughly
% the size of the estimated image displacement.  Larger values of sigma
% will make this a better approximation.  However, too large a value
% of sigma will blur out the image details and make the matching less
% precise.  sigmaBlur should be a compromise between these two effects.

rhoRewarp = 0.5;  % (default) 0.5
% Determines when iterative rewarping steps are to be done (if any, see
% mxRewarpIts above).  Let dv be the current update to the motion
% estimate.  A rewarping step is called for whenever abs(dv) has
% a component larger than (rhoRewarp*sigmaBlur).  For dv below
% this threshold, the linear brightness constancy results are
% sufficiently accurate (due to the smoothing provided by the
% gaussian filter with std dev sigmaBlur).

sigmaGrey = 4; % (default) 4
% Estimated std dev of noise in image grey levels.  This
% is used to decide on an appropriate minimum amplitude for the image
% gradients.  The idea is we do not want to be using gradients dominated
% by the noise.

vm = 2*sigmaBlur; % (default) 2* sigmaBlur
% Range of velocities to plot in the motion constraint
% plots, i.e. vx and vy are plotted in the range +-vm.

figScl = 1; % (default) 1
% For display only, blow up image figure by a factor of figScl (1 or 2).

%%%  Tracker Setup

% Image Sequence Selection  

if (chooseSeq == 1)
  %% FLEETFACE sequence
  imRoot = [matlabVisRoot '/images/seq/fleetface1/']; 
  frameStart = 250;
else
  %% DUDEKFACE sequence
  imRoot = [matlabVisRoot '/images/seq/dudekface1/'];
  frameStart = 200;
end
fRoot = 'frameDec';  % Prefix of image filename.

% Preview Image Sequence
if SEQ_PREVIEW
  % Loop over the first 51 image filenames in fn.
  figure(1); close 1;
  for frameNum = frameStart + (0:50)

    % Read the frameNum^th image
    fname = [imRoot fRoot num2strPad(frameNum, 4) '.pgm'];
    if ~exist(fname, 'file')
      error(['Cannot find file' fname]);
    end
    im = pgmRead([imRoot fRoot num2strPad(frameNum, 4) '.pgm']);
    
    % Display the frameNum^th image
    h = figure(1); showIm(im);
    resizeImageFig(h, size(im), figScl);
    pause(0.1);
    
  end
end

% Initialize Blur and Derivative Filters 

% Filter length as a function of the sigma.
gBlurSize = 2 * round(2.5 * sigmaBlur) + 1;  

% Construct filters for blurring and differentiating images.
x = (1:gBlurSize) - round((gBlurSize+1)/2);  
gFilt = exp(- x .* x / (2.0*sigmaBlur*sigmaBlur));
gFilt = gFilt/ sum(gFilt(:));         % Blur filter
gxFilt = (-x/sigmaBlur^2) .* gFilt;   % Derivative filter
filtLen = max([length(gFilt), length(gxFilt)]);

% Estimate the minimum gradient to be used in motion constraints.
% Std dev of gradient due to noise alone
gradStdDev = sigmaGrey * norm(gxFilt)*norm(gFilt); 
% minimum square of gradient norm is (2.0 std devs)^2
gTol2 = (2.0 * gradStdDev)^2; 

%%% Initialize Tracking Box

% We will track a region within a bounding box.  Here we initialize
% this bounding box.

% Read initial image frame
frameNum  = frameStart;
fname=[imRoot fRoot num2strPad(frameNum, 4) '.pgm'];
im = pgmRead([imRoot fRoot num2strPad(frameNum, 4) '.pgm']);

% Display initial image
h = figure(1); clf; showIm(im);
resizeImageFig(h, size(im), figScl);
pause(0.1);

% Select initial region
if ~USE_MOUSE 
  %% Use the regions we have selected for you
  if findstr('fleetface1',imRoot) && frameStart == 250
    mouse_pts = [83   122; 129   178];
  elseif findstr('dudekface1',imRoot) && frameStart == 0
    mouse_pts = [ 63    53; 125   133];
  elseif findstr('dudekface1',imRoot) && frameStart == 200
    mouse_pts = [    195    55; 257   135];
  else
    USE_MOUSE=TRUE;
  end
end
  
% Mouse in initial region
if USE_MOUSE
  fprintf(2, 'Select top-left and bottom-right corners of tracking region\n');
  mouse_pts = ginput(2);
end

% Display initial region
bbox = [floor(mouse_pts(1,:)); ceil(mouse_pts(2,:))];
if any(bbox(1,:)> bbox(2,:))
  display(bbox);
  error('Points must in top-left to bottom-right order');
end
figure(1); hold on;
poly = [bbox(1,:) ; bbox(1,1) bbox(2,2); bbox(2,:); ...
        bbox(2,1) bbox(1,2);  bbox(1,:)]';
drawPolys(figure(1), poly);

%%% RESTART TRACKING HERE 

% If you wish to restart the tracking, with a box you've moused in,
% then you can just paste in the code below here.  
bbox = [floor(mouse_pts(1,:)); ceil(mouse_pts(2,:))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Initialization for tracking  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set image cropping region in  first frame.
% Enlarge bounding box to have borders at least sigmaBlur wide.
% Also, the convolution code requires that the cropBox be at least
% as large as the filter width.
cropBox0 = setCropBox(bbox, sigmaBlur, filtLen, size(im));

% Grab some memory for least squares.
A = zeros(2,2); b = zeros(2,1);

% Initial motion estimate.
v = [0 0];

% Prepare to record tracking results.
if exist('trackResults', 'var')
  clear trackResults;
end
if exist('trackInit', 'var')
  clear trackInit;
end
trackInit.imRoot = imRoot;
trackInit.frameStart = frameStart;
trackInit.bbShape = bbox - repmat(bbox(1,:),2,1);
trackInit.v0 = v;
trackInit.bbCorner = bbox(1,:);
trackInit.sigmaBlur = sigmaBlur;
trackInit.mxRewarpIts = mxRewarpIts;

%%%%%%%%%%%%%%%%%%%%%%%%%% Tracking Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over the image filenames in fn.
maxFrames = 151; % 5000; (default)
for frameNum = (frameStart-1) + (1:maxFrames);

  % Read the frameNum^th image
  fname = [imRoot fRoot num2strPad(frameNum, 4) '.pgm'];
  if ~exist(fname, 'file')
    fprintf(2, ['Cannot find file ' fname '\n']);
    break;
  end

  % Plot the frameNum^th image
  im1 = pgmRead([imRoot fRoot num2strPad(frameNum, 4) '.pgm']);
  figure(1); clf; showIm(im1);
  if figScl ~= 1
    resizeImageFig(figure(1), size(im), figScl);
  end
  hold on;
  if frameNum == frameStart
    poly = [bbox(1,:) ; bbox(1,1) bbox(2,2); bbox(2,:); ...
            bbox(2,1) bbox(1,2);  bbox(1,:)]';
    drawPolys(figure(1), poly);
    pause(0.1);
  end

  % Predict the motion params
  if ZERO_PREDICT
    v = [0 0];  
  end
  % Otherwise v is given by the previous pair of frames.
  
  % Use the motion prediction for prewarping the current frame
  % We only use integer shift prewarps.
  vShift = round(v);
  % Residual from prewarp.
  dv = v-vShift;
  
  % Estimate solution of linearized brightness constancy.
  % Iterative rewarping is used if mxRewarpIts > 1
  for rewarpIts = 1: mxRewarpIts
    
    % If the discrete prewarp that was used in the previous iteration
    % is small relative to sigmaBlur, then the linear approximation
    % in the brightness constancy constraints is probably ok.
    % So we might can break out of the iterative rewarping loop early.
    if rewarpIts > 1 
      if (max(abs(dv)) < 0.6) || (norm(dv) < rhoRewarp*sigmaBlur)  
        % First condition above is used to reduce chattering...i.e. hopping
        % back and forth between two neighbouring discrete shifts.
        break;
      end
    end
    
    % Compute the shift to use for the prewarp.
    vShift = round(v);
    
    % Crop image at current frame (im1).  Note the cropBox has
    % been shifted by vShift.  This is the prewarp on the
    % first iteration, and the rewarps later.
    im1C = cropImage(im1, cropBox0 + repmat(vShift,2,1));
    
    % Blur the cropped image.
    im1CB = rconv2sep(im1C, gFilt, gFilt);
    
    % If we have seen at least two frames...
    if frameNum > frameStart
      
      % Compute the bounding box in the cropped image 
      % for the previous frame (im0).
      bboxC = bbox - repmat(cropBox0(1,:), 2,1);
      % Compute the logical image for the inside of the crop box.
      [xImC yImC] = meshgrid(1:size(im0CB,2), 1:size(im0CB,1));
      idImC = (xImC >= bboxC(1,1)) & (xImC <= bboxC(2,1));
      idImC = idImC & (yImC >= bboxC(1,2)) & (yImC <= bboxC(2,2));
      
      % Compute the linearized brightness constancy constraints
      C = linBCConstraints(im1CB, im0CB, gradIm0C, idImC, xImC, yImC, ...
                           gTol2);
      
      % Compute the least squares estimate of the update dv to the
      % translational velocity.
      if size(C,2)>0
        A(1,1) = sum(C(1,:).*C(1,:),2);
        A(2,1) = sum(C(2,:).*C(1,:),2);
        A(1,2) = A(2,1);
        A(2,2) = sum(C(2,:).*C(2,:));
        b(1) =  sum(C(1,:).*C(3,:));
        b(2) =  sum(C(2,:) .* C(3,:));
        
        iA = inv(A);
        dv = - (iA * b)';
      else
        dv = [0 0];
      end

    % Record estimates.
    if rewarpIts == 1
      frm = frameNum - frameStart;
      trackResults{frm}.vShift = vShift;
      trackResults{frm}.dv = dv;
    else
      trackResults{frm}.vShift = [trackResults{frm}.vShift; vShift];
      trackResults{frm}.dv = [trackResults{frm}.dv; dv];
    end

    % Update the velocity estimate.
    v = vShift + dv;
    fprintf('frame %d:   dv: %f %f   v: %f %f\n', frameNum, dv, v);

      % Plot the linearized brightness constancy constraints
      if PLOT_BCC 
        if rewarpIts == 1
          figure(2); clf;
        end
        kPlot = min(rewarpIts, 3);
        figure(2); subplot(1,min(3,mxRewarpIts),kPlot);
        plotLinBCC(C, vm, nPlotBCC);
        xlabel('dVx (pixels/frame)');
        ylabel('dVy (pixels/frame)');
        title(sprintf('LBCC vS: %d %d v: %4.1f %4.1f', vShift, v));
        hold on;
        plot(dv(1), dv(2), '*r', 'MarkerSize', 12);
      end

    end % if frameNum > frameStart
  end % rewarping iterations
  
  if frameNum > frameStart
    % Move bounding box  
    bbox(1:2,1) = bbox(1:2,1) + v(1);
    bbox(1:2,2) = bbox(1:2,2) + v(2);

    % Record new box position.
    frm = frameNum - frameStart;
    trackResults{frm}.bbCorner = bbox(1,:);

    % Plot bounding box
    figure(1); hold on;
    poly = [bbox(1,:) ; bbox(1,1) bbox(2,2); bbox(2,:); ...
            bbox(2,1) bbox(1,2);  bbox(1,:)]';
    drawPolys(figure(1), poly);
    
  end
  
  % Prepare for reading the next frame.
  im0 = im1; 
  cropBox0 = setCropBox(bbox, sigmaBlur, filtLen, size(im1));
  im0C = cropImage(im0, cropBox0);
  im0CB = rconv2sep(im0C, gFilt, gFilt);
  gradIm0C = zeros([size(im0CB) 2]);
  gradIm0C(:,:,1) = rconv2sep(im0C, gxFilt, gFilt);
  gradIm0C(:,:,2) = rconv2sep(im0C, gFilt, gxFilt);
  
  % Encourage Matlab to draw the figures before continuing.
  pause(0.01);
end

%%% Examine Results

% If you wish to save the results (change the filenames)
% save 'dudekfaceDefaultRes' trackInit trackResults;
% save 'fleetfaceDefaultRes' trackInit trackResults;

%load('dudekfaceDefaultRes');
%load('fleetfaceDefaultRes');
%
% Unpack trackResults cell array to get:
%   v1, the motion estimate from the linearized brightness constancy
%       constraints (LBCC) 
% and
%   v,  the motion estimate from iterative rewarping.
clear v v1 dv dv1;
fS = trackInit.frameStart;
mxRI = trackInit.mxRewarpIts;
sigB = trackInit.sigmaBlur;
for frm = 1:length(trackResults)
  v1(:, frm) = (trackResults{frm}.vShift(1,:) + ...
                trackResults{frm}.dv(1,:))';
  v(:,frm) = (trackResults{frm}.vShift(end,:) + ...
              trackResults{frm}.dv(end,:))';
  dv1(:, frm) = (trackResults{frm}.dv(1,:))';
  dv(:,frm) = v(:,frm) - (trackResults{frm}.vShift(1,:))';
end

% Plot the norm of the velocities estimated using LBCC, and those
% estimated using iterative rewarping. The difference roughly indicates
% the error in using the linear brightness constancy constraints.
nFrames = length(trackResults);
figure(3); clf;
nrmV1 = sqrt(sum(v1.^2,1));
plot(fS+(1:nFrames), nrmV1, '-og'); 
if mxRI == 1
  title('Velocity Estimates: LinBCC');
elseif  mxRI > 1
  figure(3); 
  hold on;
  nrmV1 = sqrt(sum(v1.^2,1));
  nrmV = sqrt(sum(v.^2,1));
  nrmDiff = sqrt(sum((v-v1).^2,1));
  plot(fS + (1:nFrames), nrmV1, '-og'); hold on;
  plot(fS + (1:nFrames), nrmV, '-ob'); 
  plot(fS + (1:nFrames), nrmDiff, '-r'); 
  title('Velocity Estimates: LinBCC(g), ReWarp(b), Diff(r)');
end
% Crop axes
ax = axis;
axis([fS fS+nFrames 0 min(ax(4), 8*sigB)]);
xlabel('Frame Number');
ylabel('Norm of Velocity Estimate');
% You may need to crop the plot if there are outliers in
% the motion estimates.

% Plot the norm of the velocity update from the initial rewarp
% prediction (i.e. vShift(1,:), this will be [0 0] if ZERO_PREDICT
% is true) and the difference in the velocity update using LBCC and
% iterative rewarping.  This difference approximates the error
% in using the LBCC.
if mxRI > 1
  figure(4); clf;
  plot(sqrt(sum(dv.^2, 1)), sqrt(sum((dv-dv1).^2,1)), '.b');
  ax = axis;
  axis([0 min(ax(2),8*sigB) 0 min(ax(4), 8*sigB)]);
  title('Estimated Error in LBCC');
  xlabel('Norm of Velocity Update');
  ylabel('Norm of Difference in Updates');
end
% You may need to crop the plot if there are outliers in
% the motion estimates.
%
% Due to errors in the linear approximation, the estimated error in the LBCC
% should increase (roughly quadratically) when the norm of the velocity update
% increases beyond about sigmaBlur.  (This is clearest when motion prediction
% is NOT used (i.e. ZER0_PREDICT == TRUE).  Otherwise the size of
% the velocity updates be quite small.)

%%% Other Tracking Examples

if FALSE  % Comment out the following code. But feel free to cut and paste.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% TRACKING SMALL REGIONS %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Try just tracking a small (eg. 16x16) region (eg. the tear duct on
  % one of Fleet's eyes, or the end of his nose).  To do this, set
  USE_MOUSE = TRUE; PLOT_BCC = FALSE; ZERO_PREDICT = FALSE;
  chooseSeq=1; figScl=2;
  % and mouse in the desired region in the first frame.
  % Re-execute beginning at "Image Sequence Selection" above.

  % Alternatively, try tracking a key on the keyboard in
  % the dudekface sequence (eg the left-arrow key near the bottom
  % right of the keyboard remains unoccluded). Or try tracking
  % a corner of the light switch.  Use:
  USE_MOUSE = TRUE; PLOT_BCC = FALSE; ZERO_PREDICT = FALSE; 
  chooseSeq=2; figScl=2;

  % Also try regions with less variation in the grey-levels, such
  % as a letter on the T-shirt. Or the logo/name on the top right/left
  % of the monitor (again in the dudekface sequence).
  
  % The smaller the tracking region, the more signficant the
  % accumulated errors are going to be relative to the size of
  % the box.  Currently our bounding box has size (in x and y directions):
  fprintf(2, 'Tracking box of size %3.1f by %3.1f\n', ...
          trackInit.bbShape(2,:));
  % Also, for a smaller bounding box, the errors in the
  % motion estimates can be be expected to be larger since fewer
  % constraints can be used.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TRACKING: REGIONS MUCH LARGER THAN OBJECT %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Select a region that contains the face but is about three times
  % as high and wide as the face in the first frame.  Use
  USE_MOUSE = TRUE; PLOT_BCC = TRUE; ZERO_PREDICT = FALSE;
  chooseSeq=1; figScl=1;
  % With such a sloppy region and least squares, the estimated motion
  % is influenced by the background motion.
  % Ditto for the dudekface sequence
  USE_MOUSE = TRUE; PLOT_BCC = TRUE; ZERO_PREDICT = FALSE;
  chooseSeq=2; figScl=1;
  % The estimated motion is a rough compromise between the motion of
  % the background and the motion of the head.
 
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FAILURE OF BRIGHTNESS CONSTANCY %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Try tracking the whole computer screen in the dudekface sequence.
  USE_MOUSE = TRUE; PLOT_BCC = TRUE; ZERO_PREDICT = FALSE;
  chooseSeq=2; figScl=1;
  % There are a lot of outliers in the BCC plots.
  % Also, try tracking:
  %   - just the top-left (as you look at it) quadrant of the
  %     computer screen.  Or a region half again as big, or 1/4 as big
  %   - a small (eg 16x16) region within the screen, away from the
  %     borders.
  % You should see the tracking region drift more and more as
  % the tracking region is made smaller and as it is moved
  % away from the stable screen borders. Also, the track is
  % better when there are stable (i.e. moving with the
  % screen), high contrast edges in the tracking box.  The
  % failure of brightness constancy is less severe near these
  % high contrast edges.
 
end

%%% Discussion 

% Appearance Memory.  The tracker has a very short memory (1 frame) for
% what it is tracking.  In particular, this tracker is simply tracking
% whatever image structure ends up in the bounding box (bbox) at the
% current frame.  This structure can rapidly change over time (see
% the end of the fleetface sequence), and the tracker will simply
% forget what it was tracking before, and start tracking the current
% contents of the bbox.

% Adaptation. Due to the short memory, the tracker can rapidly adapt
% to changes in the object appearance.

% Drift.  Due to the short memory, and sccumulating errors, the tracker
% can drift off of the intended target over time.  For example
% suppose the tracking error is a normal random variable with a 
% standard deviation of about sigErr pixels per frame.  How
% much drift would you expect to accumulate due to this error alone,
% say for a sequence of length nFrame? 

sigErr = 0.1;
nFrame = 1000;
figure(5); clf; hold on;
for k = 1:5
  err = sigErr * randn(2, nFrame);
  accumErr = cumsum(err,2);
  plot(1:nFrame, sqrt(sum(err.^2,1)), 'r'); 
  hold on;
  plot(1:nFrame, sqrt(sum(accumErr.^2,1)), 'b');
end
title('Simulated Drifts: Tracking Error(b), Motion Error(r)');
xlabel('Frame Number');
ylabel('Error in pixels or pixels/frame');
% Can repaste the above code into the matlab window a few
% times to get an idea of the expected variation.

% The above model of the drift is a (Gaussian) random walk.  The
% accumulated error after n frames is a Gaussian random variable
% with mean zero and variance n * sigErr^2 in each of the x,y
% components.

% We can superimpose curves indicating 1 and 2 standard deviations
% of this model on the previously sampled drifts.
figure(5); hold on;
n = 1:nFrame;
plot(n, sqrt(2 *n * sigErr^2), 'g');  % The factor of 2 is to account
                                      % for the x and y components.
plot(n, 2*sqrt(2 *n * sigErr^2), 'g');
grid on;

% For example, for sigErr = 0.1 the norm of the drift after 1000 frames
% has a standard deviation of:
fprintf(2, ' %3.1f pixels\n', 0.1 * sqrt(2 * 1000) ); 
% If we only estimated the motion to the nearest pixel, then
% the after 1000 frames, we might expect accumulated errors with a
% standard deviation of 
fprintf(2, ' %3.1f pixels\n', 1.0 * sqrt(2 * 1000) ); 
% So clearly, in order to track over 1000 frames or so, we need to keep
% the standard deviation of motion errors well below 1 pixel/frame.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% End: Motion Tutorial %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
