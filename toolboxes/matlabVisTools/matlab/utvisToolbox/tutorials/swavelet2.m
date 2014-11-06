% swavelet2.m  STEP 2 of self-inverting multiscale pyramid.
% 
% Demo of self-inverting multiscale pyramid, as designed in
% CSC 2503, Lecture 4.98.
%
% To run this script, you need to have the filter kernels
% set up according to step 1 (see swavelet1.m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 0:
%%  Check preconditions...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you are running the student edition of Matlab 
% (i.e. your prompt says EDU>>) then you have a limit on the
% size of arrays.  You should put the line 
% SIZE_LIMIT = 1;
% in your startup.m file (uncommented...remove the %)
% and also execute it this line now.
  
delay = 1;
if (~(exist('lp0','var') & exist('hp0','var') ...
      & exist('lp','var') & exist('hp','var')))
 fprintf(1, 'Heh, I said run Step 1 first!\n');
end

fprintf(1, 'Filter radii %d %d\n', fS2, f0S2);
if (max([fS2 f0S2]) > 15)
 fprintf(1, 'Warning: Filter radii larger than expected.\n');
 fprintf(1, '...      Filtered images must be bigger than the filters.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2: Apply filter kernels (l0(x),h0(x)), (l(x), h(x)) in a
%%         cascade.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Input an image...
im = pgmRead('einstein.pgm');
%% This particular image has a bothersome border, remove it
sizeIm = size(im);
border = 2;
im = im(border+1:sizeIm(1)-border, border+1:sizeIm(2)-border);
clf;
showIm(im);

%% Correlate the image with filters lp0 and lp1.  (Correlation
%% is just convolution, but without reflecting the filter
%% kernels about the origin in x,y.)
%% See help corrDn for a description of how it handles
%% the borders of the image. (There are several options.)
border_type = 'reflect1';
im_lp0 = corrDn(im, lp0, border_type);
im_hp0 = corrDn(im, hp0, border_type);

% Peek at results
showIm(im_lp0); pause(delay);
showIm(im_hp0); pause(delay);

% Check FT
qw = abs(fftshift(fft2(im_hp0)));
showIm(log(1.0+qw)); pause(delay);

qw = abs(fftshift(fft2(im_lp0)));
showIm(log(1.0 + qw)); pause(delay);

% Try reconstructing from just these two filtered images (lp0,hp0)
% by convolving with these filters (instead of correlating) and add
% the results...

im_recon = upConv(im_lp0, lp0, border_type);
im_recon = im_recon + upConv(im_hp0, hp0, border_type);
reconErr = im_recon-im;

% Check out the reconstruction.
showIm(im_recon); pause(delay);
imStats(im, im_recon);

% Show the residual image (BEWARE it has been rescaled for display).
showIm(reconErr); pause(delay);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check out lp and hp cascade.
% This is written in a very simple way to ensure
% that it runs on a variety of older versions of Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_lp0 = corrDn(im, lp0, border_type);
im_hp0 = corrDn(im, hp0, border_type);
sz = size(im);

im_bp1 = corrDn(im_lp0, hp, border_type);
im_lp = corrDn(im_lp0, lp, border_type, [2 2]);
sz = [sz ; size(im_lp)];

im_bp2 = corrDn(im_lp, hp, border_type);
im_lp = corrDn(im_lp, lp, border_type, [2 2]);
sz = [sz ; size(im_lp)];

im_bp3 = corrDn(im_lp, hp, border_type);
im_lp = corrDn(im_lp, lp, border_type, [2 2]);
sz = [sz ; size(im_lp)];

%% Look at intermediate results.
showIm(im_hp0); pause(delay);
showIm(im_bp1); pause(delay);
showIm(im_bp2); pause(delay);
showIm(im_bp3); pause(delay);
showIm(im_lp); pause(delay);

%% Check out fft of filtered images
nZero = floor((sz(1,:)+1)/2)+1;
% Choose one of the following ffts and displays:
 qw = abs(fftshift(fft2(im_hp0)));
 %% Display log of amplitude spectrum (for values above 1000 or so)
 qw = log(1000)-log(qw+1000);
 showIm( qw ); pause(delay);
 qw = abs(fftshift(fft2(im_bp1)));
 qw = log(1000)-log(qw+1000);
 showIm( qw ); pause(delay);
 qw = embed(abs(fftshift(fft2(im_bp2))), sz(1,:));
 qw = log(1000)-log(qw+1000);
 showIm( qw ); pause(delay);
 qw = embed(abs(fftshift(fft2(im_bp3))), sz(1,:));
 qw = log(1000)-log(qw+1000);
 showIm( qw ); pause(delay);
 qw = embed(abs(fftshift(fft2(im_lp))), sz(1,:)); qw(nZero, nZero) = 0.0;
 qw = log(1000)-log(qw+1000);
 showIm( qw ); pause(delay);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct image by upsampling and convolving low pass versions
% of the image with lp and bandpass versions with hp. These are
% the same filters used above for the high-, band-, and low-pass
% decomposition.  Hence the filter bank is said to be 'self-inverting'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_recon = upConv(4 * im_lp, lp, border_type, [2 2], [1 1], sz(3,:));
im_recon = im_recon + upConv(im_bp3, hp, border_type);
% The factor of 4 above is necessary.  This comes from the change
% in the effective area of each discrete pixel.  When pixels are subsampled
% by 2 in each direction, x and y, they cover 4 times the area.

im_recon = upConv(4 * im_recon, lp, border_type, [2 2], [1 1], sz(2,:));
im_recon = im_recon + upConv(im_bp2, hp, border_type, [1 1]);

im_recon = upConv(4 * im_recon, lp, border_type, [2 2], [1 1], sz(1,:));
im_recon = im_recon + upConv(im_bp1, hp, border_type);

im_recon = upConv(im_recon, lp0, border_type);
im_recon = im_recon + upConv(im_hp0, hp0, border_type);

%% Check out error in reconstructed image: im_recon

borderWidth = 0;
sizeIm = size(im);
reconErr = im_recon-im;
clf;
showIm(reconErr(1+borderWidth:sizeIm(1)-borderWidth, ...
                1+borderWidth:sizeIm(2)-borderWidth)); pause(delay);

range(reconErr)
          
imStats(im, im_recon)

showIm(im_recon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Completed self-inverting pyramid demo.
%%
%% Things to do:
%% a) Accuracy.
%%    Try changing the filter sizes (fS2, f0S2), and check out effect on
%%    reconstruction.
%%    We saw in part 1 that larger filter sizes allowed for a more accurate
%%    partitioning of the power spectrum by lp0, hp0, and lp, hp.  This
%%    should translate into a more accurate reconstructed image im_recon.
%%    For a fixed filter size it is likely that more accurate
%%    pairs of filters (lp0, hp0) and (lp, hp) can be found.  Here
%%    we simply truncated the inverse Fourier transform of the designed
%%    partitioning.
%% b) Interpolation. Note that the low frequency information can be 
%%    upsampled using a similar inverse procedure, but with the information
%%    for finer scales simply set to be zero.  Try this for im_lp,
%%    arriving at an image of the same size as the original im, but
%%    with only the low-pass part of the Fourier spectrum.  Similarly,
%%    interpolate the information in im_lp and im_bp3 up to the original
%%    pixel resolution.
%% c) An approximate inverse can also be useful, especially when you
%%    consider that small filter sizes can then be used.  Try expanding
%%    an image I iteratively, as in
%%    b0 = 0
%%    e0 = I
%%    for n>=0
%%        e_{n+1} = e_n - B b_n
%%        stop when e_n sufficiently small
%%        b_{n+1} = B^T e_{n+1}
%%    end
%%    Then set
%%       a = b_0 + b_1 + b_2 ... + b_N
%%    so
%%       Ba = B b_0 + ... + B b_N
%%          = Sum( e_n - e_n+1) = e_0 - e_N = I - e_N
%%    This loop will converge so long as B B^T is close enough
%%    to the identity (i.e. the eigenvalues of BB^T -I are within
%%    the unit circle). 
%%    It provides an accurate set of expansion coefficients, a,
%%    such that Ba is a good approximation for I (say within 1/2 a
%%    grey-level).
%% d) A similar partitioning of the power spectrum can be done in
%%    orientation.  Use half cosine windows in the angle of the
%%    frequency (k_x, k_y).  Ensure that the power of these cosine
%%    windows sums to one.  Apply this decomposition to each of
%%    the scales of the current self-inverting multiscale pyramid.
%%    (See also the papers by Simoncelli on the reading list.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
