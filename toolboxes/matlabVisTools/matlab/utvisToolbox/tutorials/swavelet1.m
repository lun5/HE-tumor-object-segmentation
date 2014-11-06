% swavelet1.m:
%  Script for STEP 1 of self-inverting multiscale pyramid demo.
%
% Demo of self-inverting multiscale pyramid, as
% designed in CSC 2503, Lecture 4.98.  (Do not run this file,
% but rather open it with an editor, and paste commands into
% matlab.

%%% STEP 0: The scale-wavelets use some M-files in the
%%  tutorials/localUtils directory.  Add this to your path.

global matlabVisRoot

% We need to ensure the path is set for the iseToolbox.
if isempty(matlabVisRoot)
  dir = pwd;
  cd ~jepson/pub/matlab   % CHANGE THIS
  startup;
  cd(dir);
end
addpath([matlabVisRoot '/utvisToolbox/tutorials/localUtils']);
which embed  % should find it in tutorials/localUtils

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: Build filter kernels (l(x),h(x)), (l0(x), h0(x))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% We will build the 2D Fourier transform of the required
% filter pairs (l,h), and (l0,h0) with N samples in each dimension
N = 128;

%%%%%%%%%%%
%% STEP 1a: Build the 2D Fourier transforms of the filters
%%%%%%%%%%%
% Build a 1D radial cross-section of the transform of each filter
ffilt = fcos_wave(N);
nZero = floor(N/2)+1;

% Take a peek at them
hold off;
plot(ffilt'); axis([1 N 0 1.1]);

% Build 2D versions using interpolation... 
ramp = mkR([N N], 1, [nZero, nZero]);
flp0 = pointOp(ramp,[ffilt(3,[nZero:-1:1]) 0.0], 0, 1, 0);
fhp0 = pointOp(ramp,[ffilt(4,[nZero:-1:1]) 1.0], 0, 1, 0);

flp = pointOp(ramp,[ffilt(1,[nZero:-1:1]) 0.0], 0, 1, 0);
fhp = pointOp(ramp,[ffilt(2,[nZero:-1:1]) 1.0], 0, 1, 0);

clear ramp;

% Take a peek at them (one at a time)
delay = 1;
showIm(flp); pause(delay);
showIm(fhp); pause(delay);
showIm(flp0); pause(delay);
showIm(fhp0); pause(delay);

% Plot the diagonal of one of these...
plot(diag(fhp0)); axis([1 N 0 1.1]);
% notice how it has been extrapolated to be constant at either end

%%%%%%%%%%%
%% STEP 1b: Take the inverse FT, and crop the result to form
%%          a fS x fs filter kernel, where fS = 2*fS2 + 1
%%%%%%%%%%%

% Set the radius of the desired filter kernels
% Larger values here give more accurate self-inversion (that
% is, the cropped filters are better approximations to the designed
% FTs).  However, the filtering will take longer (time proportional
% to the fS^2).

f0S2 = 10;  % filter radius of lp0, hp0 pair
fS2 = 10;   % filter radius of lp, hp pair  (need not be the same)
% 10 here is good for starters... but try values 5 10 15 later

% Build truncated (lp0, hp0) pair:
lp0 = real(fftshift(ifft2(fftshift(flp0))));
lp0 = lp0(nZero-f0S2: nZero+f0S2, nZero-f0S2:nZero+f0S2);
hp0 = real(fftshift(ifft2(fftshift(fhp0))));
hp0 = hp0(nZero-f0S2:nZero+f0S2, nZero-f0S2:nZero+f0S2);
% Make sure high-pass filter sums to zero.
hp0(f0S2+1,f0S2+1) = hp0(f0S2+1, f0S2+1) - sum(sum(hp0));
% Rescale cropped filters to fit energy partition as well as possible
[lp0 hp0] = adj_partition(lp0, hp0, N);

% Build truncated (lp, hp) pair
lp = real(fftshift(ifft2(fftshift(flp))));
lp = lp(nZero-fS2: nZero+fS2, nZero-fS2:nZero+fS2);
hp = real(fftshift(ifft2(fftshift(fhp))));
hp = hp(nZero-fS2:nZero+fS2, nZero-fS2:nZero+fS2);
% Make sure high-pass filter sums to zero.
hp(fS2+1,fS2+1) = hp(fS2+1, fS2+1) - sum(sum(hp));
% Rescale cropped filters to fit energy partition as well as possible
[lp hp] = adj_partition(lp, hp, N);

% Take a peek at the filter kernels
clf;
subplot(1,2, 1), showIm(lp0), subplot(1,2, 2), mesh(lp0), pause(delay);
subplot(1,2, 1), showIm(hp0), subplot(1,2, 2), mesh(hp0), pause(delay);
subplot(1,2, 1), showIm(lp), subplot(1,2, 2), mesh(lp), pause(delay);
subplot(1,2, 1), showIm(hp), subplot(1,2, 2), mesh(hp), pause(delay);

%%Check energy partition of cropped filter masks
alp0 = abs(fftshift(fft2(lp0,N,N))).^2;
ahp0 = abs(fftshift(fft2(hp0,N,N))).^2;
alp = abs(fftshift(fft2(lp,N,N))).^2;
ahp = abs(fftshift(fft2(hp,N,N))).^2;

% Display power spectra for cropped filters
clf;
showIm(alp0); pause(delay);
showIm(ahp0); pause(delay);
showIm(alp); pause(delay);
showIm(ahp); pause(delay);

% Pairs of these should add to 1...roughly.
% Compute and display this sum of squares.
% BEWARE the sum has been rescaled for display purposes, check
% out the range.
qw = alp0 + ahp0;
range(qw)
showIm(qw); pause(delay);

qw = alp + ahp;
range(qw)
showIm(qw); pause(delay);


% More informative is the plot of the diagonal or a row through the middle:
% Choose one of:
 fit = alp0; model = flp0.^2;  % Remember to square the model amp. spectrum.
 fit = ahp0; model = fhp0.^2;
 fit = alp0 + ahp0; model = ones(N);  % Or model = flp0.^2 + fhp0.^2;
 fit = alp + ahp; model = ones(N);  % Or model = flp.^2 + fhp.^2;
 fit = alp; model = flp.^2;  
 fit = ahp; model = fhp.^2;
% And choose one of:
 plot(diag(fit),'r'); axis([1 N 0 1.1]);
 hold on; plot(diag(model),'g'); axis([1 N 0 1.1]);
% .. or ..
 plot(fit(nZero,:),'r'); axis([1 N 0 1.1]);
 hold on; plot(model(nZero,:),'g'); axis([1 N 0 1.1]);
hold off;

%%%%%%%%%%% PLAY %%%%%%%%%%%%%%%%%%%
%    Try the maximum possible size: 
%     fS2 = N-nZero;
%     f0S2 = N - nZero;
%    Redo the above, going back to step 1b (just after setting fS2). 
%    You should get a clean fit.  This shows the error is
%    due to truncating the filter kernel to be fS x fS.
%    Try intermediate sizes, eg: fS2 = 30;
%    Make sure you reset fS2 to a more reasonable
%    value before continuing on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear fit model qw;
clear ahp alp ahp0 alp0 ffilt fhp fhp0 flp flp0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finished Step 1:  
%%   Open swavelet2.m to try these filters out.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
