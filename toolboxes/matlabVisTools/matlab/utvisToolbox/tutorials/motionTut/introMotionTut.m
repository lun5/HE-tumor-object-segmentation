%% NOTE: make sure you first go to the existing motion tutorial directory
% cd to whereever you need to be.

clear
close all
global matlabVisRoot  % Root directory for ise and utvis toolboxes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Let's first grab two frames of an image sequence

% this is where I cam grabbing the images from.  You may
% want to change this and keep the images locally.
imRoot = [matlabVisRoot '/images/seq/fleetface1/'];
fRoot = 'frameDec';  % Prefix of image filename.
frameNum = 250;
im0 = pgmRead([imRoot fRoot num2strPad(frameNum, 4) '.pgm']);
im1 = pgmRead([imRoot fRoot num2strPad(frameNum+1, 4) '.pgm']);
imSz = size(im0);

%% Let's display the images to get a sense of the motion.
figure(1);
for n = 1:8
   showIm(im0); pause(1/10);
   showIm(im1); pause(1/10);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Now, let's choose a rectangular region within which we'll estimate 
%   the motion, and display it on top of the image in the first frame
bbox = [83   122; 129   178];
figure(1); hold on;
showIm(im0);
resizeImageFig(figure(1), imSz, 2);
% showIm(im1);  ????????
poly = [bbox(1,:) ; bbox(1,1) bbox(2,2); bbox(2,:); ...
        bbox(2,1) bbox(1,2);  bbox(1,:)]';
drawPolys(figure(1), poly);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  We'll need gradient filters: We've already defined these offline:

sigmaBlur = 2; 			% Gaussian Sigma for all filters.
sigmaGrey = 4;			% Estimated std dev of noise (in grey levels)  
%  Filter length as a function of the sigma.
gBlurSize = 2 * round(2.5 * sigmaBlur) + 1; 
%  Construct filters for blurring and differentiating images.
x = [1:gBlurSize] - round((gBlurSize+1)/2);  
gFilt = exp(- x .* x / (2.0*sigmaBlur*sigmaBlur));
gFilt = gFilt / sum(gFilt(:));         			% Blur filter
gxFilt = (-x/sigmaBlur^2) .* gFilt;   			% Derivative filter
filtLen = max([length(gFilt), length(gxFilt)]);

% Estimate min gradient to be used in motion constraints.
% Std dev of gradient due to noise alone
gradStdDev = sigmaGrey * norm(gxFilt)*norm(gFilt); 
% minimum square of gradient norm is (2.0 std devs)^2
gTol2 = (2.0 * gradStdDev)^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Since we're only working in a small region, let's crop the image
%  (but of course it will have to be 'enlarged' to allow for the 
%  fact that we are filtering the image).
cropBox0 = setCropBox(bbox, sigmaBlur, filtLen, imSz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Now let's break down the motion estimation steps.

% First let's crop the two images, since we only need to 
% process the region of interest we specified above.  
im0C = cropImage(im0, cropBox0);
im1C = cropImage(im1, cropBox0);

% Now blur the images and take its partial derviatives.
% It is important that we use a set of consistent filters.... i.e.,
%  the derivative filters are really the derivative of a Gaussian, and
%  we need to make sure that the same Gaussian is applied to the images
%  before the temporal difference is taken.  Otherwise the brightness
%  constancy assumption may not hold (i.e. we need to take spatial derivatives
%  and temporal differences of the same image signal).
im0CB = rconv2sep(im0C, gFilt, gFilt);
im1CB = rconv2sep(im1C, gFilt, gFilt);
% We only need to linearize one of the two images
gradIm0C = zeros([size(im0CB) 2]);
gradIm0C(:,:,1) = rconv2sep(im0C, gxFilt, gFilt);
gradIm0C(:,:,2) = rconv2sep(im0C, gFilt, gxFilt);

% Let's plot the region of interest in the two images, along with
% the filter outputs.
figure(1); clf;
for n = 1:8
   showIm(im0C, 'auto', 4); pause(1/10);
   showIm(im1C, 'auto', 4); pause(1/10);
end
figure(2)
showIm(im0C + i*im1C, 'auto', 4, 'Cropped images from frame 0 and frame 1');
figure(3)
showIm(im0CB + i*im1CB, 'auto', 4, 'Blurred (low-pass) cropped images');       
figure(4)
showIm(gradIm0C(:,:,1) + i*gradIm0C(:,:,2) , 'auto', 4, 'Horizontal (left) and vertical (right) image derivatives');


% Compute the bounding box in the cropped image for the previous frame (im0).
bboxC = bbox - repmat(cropBox0(1,:), 2,1);
% Compute the logical image for the inside of the crop box.
[xImC yImC] = meshgrid(1:size(im0CB,2), 1:size(im0CB,1));
idImC = (xImC >= bboxC(1,1)) & (xImC <= bboxC(2,1));
idImC = idImC & (yImC >= bboxC(1,2)) & (yImC <= bboxC(2,2));

%% Extract the linearized brightness constancy constraints and place 
% them in a large matrix, with one constraint per row consisting 
% of the image location and the image gradient.
C = linBCConstraints(im1CB, im0CB, gradIm0C, idImC, xImC, yImC, gTol2);
      
% Now form the normal equation matrix and vector (i.e., for Av = b)
A = zeros(2,2); b = zeros(2,1);
A(1,1) = sum(C(1,:).*C(1,:),2);
A(2,1) = sum(C(2,:).*C(1,:),2);
A(1,2) = A(2,1);
A(2,2) = sum(C(2,:).*C(2,:));
b(1) =  sum(C(1,:).*C(3,:));
b(2) =  sum(C(2,:) .* C(3,:));
% Finally, compute the least squares estimate of translational velocity.
iA = inv(A);
dv = - (iA * b)';

% Our initial guess was [0, 0] so our first velocity estimate is dv
v = dv;
fprintf('frame %d, iteration 1,  dv: %f %f   v: %f %f\n', frameNum, dv, v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Now let's plot a collection of the motion constraints along with the
%  estimated velocity v.

figure(5); subplot(1, 2, 1);
vm = 2*sigmaBlur; 		% Range of vel for plotting motion constraints
nPlotBCC = 200; 		% Num of motion constraints to plot.  
plotLinBCC(C, vm, nPlotBCC);
xlabel('dVx (pixels/frame)'); ylabel('dVy (pixels/frame)');
title(sprintf('LBCC v: %4.1f %4.1f', v));
hold on;
plot(dv(1), dv(2), '*r', 'MarkerSize', 14, 'LineWidth', 2);

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Now we have a new initial guess!  So we can warp the image by a
%   new (integer) velcoity, and refine the estimate of the velocity.

vShift = round(v);	% use an integer initial guess
im1C = cropImage(im1, cropBox0 + repmat(vShift,2,1));
im1CB = rconv2sep(im1C, gFilt, gFilt);

%  Let's redisplay the images to get a sense of the motion.
figure(1); clf;
for n = 1:8
   showIm(im0C, 'auto', 4); pause(1/10);
   showIm(im1C, 'auto', 4); pause(1/10);
end
% Notice how most of the translation motion is now gone from this image pair.
% I.e., the center of the region of interest is not moving much.
% Also note that the translational model of image velocity is not a 
% very good model for the 2D motion field in this region.

%% Compute the bounding box in the cropped image for the previous frame (im0).
bboxC = bbox - repmat(cropBox0(1,:), 2,1);
% Compute the logical image for the inside of the crop box.
[xImC yImC] = meshgrid(1:size(im0CB,2), 1:size(im0CB,1));
idImC = (xImC >= bboxC(1,1)) & (xImC <= bboxC(2,1));
idImC = idImC & (yImC >= bboxC(1,2)) & (yImC <= bboxC(2,2));

C = linBCConstraints(im1CB, im0CB, gradIm0C, idImC, xImC, yImC, gTol2);
% Compute the least squares estimate of the update dv 
A(1,1) = sum(C(1,:).*C(1,:),2);
A(2,1) = sum(C(2,:).*C(1,:),2);
A(1,2) = A(2,1);
A(2,2) = sum(C(2,:).*C(2,:));
b(1) =  sum(C(1,:).*C(3,:));
b(2) =  sum(C(2,:) .* C(3,:));
iA = inv(A);
dv = - (iA * b)';

% new velocity
v = vShift + dv;
fprintf('frame %d, iteration 2,  dv: %f %f   v: %f %f\n', frameNum, dv, v);

figure(5); subplot(1, 2, 2);
vm = 2*sigmaBlur; 		% Range of vel for plotting motion constraints
nPlotBCC = 200; 		% Num of motion constraints to plot.  
plotLinBCC(C, vm, nPlotBCC);
xlabel('dVx (pixels/frame)'); ylabel('dVy (pixels/frame)');
title(sprintf('LBCC vS: %4.1f %4.1f  v: %4.1f %4.1f', vShift, v));
hold on;
plot(dv(1), dv(2), '*r', 'MarkerSize', 14, 'LineWidth', 2);

% Now the updated velocities are less than half a pixel, so we stop
% iterating the estimator with preshifts.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Stabilization:
%     OK, now let's see how well we managed to stabilize the two images
%     using the final results

im0C = cropImage(im0, cropBox0);
im1C = cropImage(im1, cropBox0);

% create vector field with velcoity v then call interp2 to shift image 1
[xi, yi] = meshgrid(1:51,1:61);
newIm1C = interp2(xi, yi, im1C, xi+v(1), yi+v(2));
idx = isnan(newIm1C(:));
if any(idx)
  newIm1C(idx) = 175;
end
  

newIm0 = [im0C, 175+zeros(size(im1C,1),20), im0C];
newIm1 = [im1C, 175+zeros(size(im1C,1),20), newIm1C];
% Let's redisplay the images to get a sense of the motion.
figure(10); clf;
for n = 1:8
   showIm(newIm0, 'auto', 4); pause(1/10);
   showIm(newIm1, 'auto', 4); pause(1/10);
end






