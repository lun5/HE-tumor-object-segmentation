%=======================================================================
% Tutorial: Introduction to Image Processing in Matlab 
%  DJF, Sept 97
%=======================================================================

% This tutorial is designed to give you a brief introduction to
% some of the basic tools and techniques for processing image in
% matlab. It also serves to outline the basic types of operations
% that are central to image processing.  The basic types of
% operations are divided into four groups for convenience:
%    1) image point operations
%    2) filtering (local operations)
%    3) statistical operations
%    4) geometric operations

% Ensure there are no name conflicts with previous work
% (First save any previous results you wish to keep.)
clear;

%-----------------------------------------------------------------------
% 1. IMAGE GENERATION AND DISPLAY

% Images are generally stored as matrices, that is, 2d arrays of
% numbers.  They may be integers, floating point numbers, or
% complex-valued numbers.  We'll make heavy use of a function
% called showIm that displays a matrix as an image.  The
% elements are called pixels (picture elements).

% There are many ways to generate interesting and useful
% synthetic images: For example, to make a bright circular disk
% on a dark background, and then display it (Cut and paste
% the following Matlab commands into the Matlab command window).
%
clf;
disc = mkDisc([64 64],16);
showIm(disc);
drawnow

% To get more information on the function makeDisc, try this:
%
help mkDisc

% To find out where it is defined, so you can go and see the code, try
%
which mkDisc

% Here, [64 64] determines the size of the image array and 16
% specifies the radius of the disk.  We can create a smaller
% image:
%
disc = mkDisc([8 8],2);
showIm(disc);
drawnow

% to show the values at each pixel, try
%
disc

% Note that we could have created a simpler matrix directly by typing:
%
disc = [0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 1 1 0 0 0;
	0 0 1 1 1 1 0 0;
	0 0 1 1 1 1 0 0;
	0 0 0 1 1 0 0 0;
	0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0];
showIm(disc);
drawnow

% But it's much easier to write functions to generate synthetic
% images.

% Another important image type is the sinusoidal image.  Here's a
% sinusoid with wavelength of 16, and an orientation of 0
% (vertical):
%
s = mkSine([64 64],16,0);
showIm(s);
drawnow

% Sometimes it's informative to see the intensity values along a
% single row in an image array... e.g., to plot the 8th row:
%
plot(s(8,:))
drawnow

% The notation s(8,:) extracts the 8th row.  s(:,3) yields the
% 3rd column.  Also sometimes we wish to convert a matrix into a
% vector by stacking the column. To do this we use s(:).

% To see how we can generate and display complex-valued images,
% let's generate a sine and cosine complex-valued image.  To do
% this we'll use more of the parameters to this function, Here's
% what help mkSine tells you about the parameters:
%
%   res = mkSine([ydim xdim], period, direction, amplitude, phase, origin)
%
s = mkSine([64 64],16,0,1.0,0.0);
c = mkSine([64 64],16,0,1.0,pi/2);
imPair = c + i*s;

% note, i is a special complex-valued number... what is i^2 ???
%
i^2

% So, best not to define this special variable as something else!

% Now, let's display the complex-valued image with the real part
% is on the left and the imaginary part is on the right:
%
showIm(imPair)
drawnow

% And let's plot a row of its values ... in the complex plane no
% less!
%
clf, plot(imPair(1,:)), axis square;
drawnow

% The function clf clears the figure window.

% If you look in the synth directory
% will see functions for generating other types of images.  Among
% these are noise generation routines.  For example:
%
im = mkFract([128,128], 0.2);
showIm(im);
drawnow

% The function truesize sets the size of the window so that each entry
% in the array is displayed at one framebuffer pixel.
% truesize;  %% I don't have the function on my laptop so it's commented out.

% A standard matlab called rand produces a matrix of uniform
% noise, and randn produces a matrix of gaussian noise:
%
noise = randn(128);
showIm(noise);
drawnow

%-----------------------------------------------------------------------
% 2. LOADING AND STORING IMAGES (from/to files)

% There are many image formats out there. A simple format is the
% pgm format (ppm for colour images).  A bunch of images should
% already available in your matlab path (take a look in the
% images subdirectory of the ise toolbox).  Here's one:
%
al = pgmRead('einstein.pgm');
showIm(al);
drawnow

% To store an image to a file, try using pgmWrite (NOTE: You'll
% have to type in your own home directory for homeDir).  To
% be safe, these commands are currently commented out:
%
% homeDir = './';
% filename = 'tmpAl.pgm';
% pgmWrite(al,[homeDir,filename]);

% pgm format only saves 8bits per pixel, with automatic (by
% default) scaling of the values to fill the range 0-255 To save
% matlab variables, including images, in full precision you can
% use the save command (Note that images will take up lots of
% disk space) (These are currently commented out).
%
% chdir(homeDir);
% save tmpAl al;

% to later load al back into memory, use the load command
% (This is currently commented out):

% load tmpAl

% to view other images, check the iseToolbox images directory.
% for example:
load aztec.mat
clf; showIm(X); colormap(map);

load parkbench.mat
clf; showIm(X); colormap(map);

ron = pgmRead('reagan.pgm');
clf; showIm(ron); colormap(gray(256));

%% Check your current storage
whos
clear X ron map

%-----------------------------------------------------------------------
% 3. SAMPLING AND QUANTIZATION

% When images are represented in this way there are two important
% parameters that govern the fidelity of the image.  One is the
% number of pixels/samples, and this is often referred to as the
% resolution of the image.  The other is the number of bits used
% to represent the intensity value at every pixel.

% To see the effect of sampling an image at different resolution,
% let's subsample Al a few times, at every 4, 8 and 16 pixels.
%
al_s4 = al(1:4:256,1:4:256);
al_s8 = al(1:8:256,1:8:256);
al_s16 = al(1:16:256,1:16:256);

% What are the sizes of these images:
%
[size(al); size(al_s4); size(al_s8); size(al_s16)]

% Let's take a look at them (Note: first make your figure window bigger):
%
close, subplot(2,2,1)
showIm(al, [0 255],0.5);
title('original')
subplot(2,2,2)
showIm(al_s4,[0 255], 2);
title('subsampled by 4')
subplot(2,2,3)
showIm(al_s8, [0 255], 4);
title('subsampled by 8')
subplot(2,2,4)
showIm(al_s16, [0 255], 8);
title('subsampled by 16')
drawnow

% The original al image has 8 bits per pixel and can therefore
% represent 256 different gray levels at each pixel.  Let's
% reduce the number of bits per pixel by dividing and rounding
% (also called quantizing):
%
al_q8 = round(al/8)*8;
al_q16 = round(al/16)*16;
al_q32 = round(al/32)*32;

% These images have 5, 4, and 3, bits per pixel. Let's take a
% look:
%
close, subplot(2,2,1)
showIm(al,[0 255], 0.5);
title('original: 8 bits/pixel')
subplot(2, 2,2)
showIm(al_q8, [0 255], 0.5);
title('5 bits/pixel')
subplot(2,2,3)
showIm(al_q16, [0 255], 0.5);
title('4 bits/pixel')
subplot(2,2,4)
showIm(al_q32, [0 255], 0.5);
title('3 bits/pixel')
drawnow

%-----------------------------------------------------------------------
% 4. IMAGE POINT OPERATIONS

% Perhaps the simplist class of operations on images are those
% that operate on pixel values directly.  So the input is an
% image and the output is an image of the same size.  The same
% operation is applied to every pixel.  The quantization that we
% just performed is an example of a point operation.

% Point operations can also operate on a pair of input images:
% for example, algebraic operations such as pixel-wise addition
% and multiplication.

% Here we add the sine and cosine images defined above:
%
s = mkSine([64 64],16,0,1.0,0.0);
c = mkSine([64 64],16,0,1.0,pi/2);
newSine = c+s;
clf
showIm(newSine);
drawnow

% It still looks like a sinusoid.  Remember your trig identities:
%      cos(kx + phi) = cos(kx) * cos(phi)  -  sin(kx) * sin(phi)
% where in our case phi=-pi/4.  Therefore, by simply taking
% linear combinations of s and c we can create sinusoidal images
% of arbitrary phase.

% But remember that cos(pi/4)=0.7071, and therefore we also
% scaled the image by 1/0.7071.  Here's another way to calculate
% the same result:
%
phi = -pi/4;
anotherSine = (c * cos(phi)  -  s * sin(phi)) / 0.7071;

% We can use a function called imStats (described below) to tell
% us about the values in the difference between two images.  If
% the images are similar then the mean and variance of the
% difference should really be close to zero.
%
imStats(anotherSine, newSine);
% or, for just the mean square error
mse(anotherSine, newSine);

% The multiplication and division operations performed above
% involved a product/quotient of a scalar and an image.  We can
% also pointwise multiply 2 images, so the product of images c
% and s at any one pixel is the product of the corresponding
% pixels at the same location in c and s.  First make a pair of
% sinusoidal grating images with different periods:
%
c1 = mkSine([64 64],64,0,1,pi/2);
c2 = mkSine([64 64],8,0,1,pi/2);

% The to multiply these images pointwise, use the operator .*
% (point multiply):
%
prod = c1 .* c2;
showIm(prod);
drawnow

% Let's look at a slice through this image
%
plot(prod(20,:));
drawnow

% If you use * instead of .* you do conventional matrix
% multiplication instead of pixel/point-wise multiplication, that
% is very different:
%
showIm(prod + i*(c1 * c2));
drawnow

% One common use of image point multiplication is windowing.
% That is, we define a signal that is nonzero in only a small
% region of the image.  Then when we multiply it with the image
% the result image to zero out structure outside the window.
%
window = mkDisc(size(al),50);
clf
showIm(al .* window);
drawnow

% There are many other image point ops, including logarithms,
% exponentiation, etc.  Take a look in the imops directory of the
% ise toolbox to see the list of them. You could implement these
% using loops, but that is way slow in matlab.  Don't do it!

% Sometimes images have a wide range of values that cannot be
% conveniently displayed on a monitor. Then one might clip the
% values.  I.e., pixels with values smaller than some threshold
% are set to the low-threshold, and values higher than another
% threshold are set to that high-threshold.
%
showIm(clip(al,100,150));
drawnow

% Thresholding is sometimes a good way to segment the image into
% useful regions, separating dark from light regions.  We'll see
% more about this later in this tutorial.

%-----------------------------------------------------------------------
% 5. IMAGE FILTERING

% Many operations involve small image regions rather than
% individual pixel values.  The relative positions of the pixels
% are therefore important.  These are called local operations.

% A simple but useful local operation is a local averaging
% (blurring) operator.  The goal is to compute the local mean
% pixel value in every small region of the image.  The local
% means form the result image. This blurring can be important to
% remove noise or irrelevant detail in an image.  A simple, but
% inefficient, way to do this is to loop over all image
% locations, creating a window centered at each point,
% multiplying it with the image.  The sum of pixel values in the
% windowed image, divided by the number of nonzero pixels in the
% window is that local mean.  E.g., the window above was centered
% at (x,y)=(180,115).  The window was 1 inside the disk. The
% function sum2 will sum the pixel values in a 2d image.  Hence
% the local average is:
%
localAve = sum2(al .* window) / sum2(window)

% A better way to do this is with convolution.  The averaging
% operator is actually called a low-pass filter.  To compute the
% average of every 7x7 region we define a filter and then apply
% it:
%
uniformFilter = (1/49) * [1 1 1 1 1 1 1; 
                          1 1 1 1 1 1 1;
			  1 1 1 1 1 1 1;
			  1 1 1 1 1 1 1;
			  1 1 1 1 1 1 1;
			  1 1 1 1 1 1 1;
			  1 1 1 1 1 1 1];
blurAl = conv2(al,uniformFilter,'same');
showIm(al + i*blurAl)
drawnow

% There are several convolution functions that you will be using.
% Conv2 a standard matlab function.  The others are in the ise
% toolbox.  The basic operation is always the same, but the ise
% toolbox functions provide additional flexibility, e.g., in how
% they handle the edges of the image.

% Let's add some noise to Al, and see what the blurring does:
%
noisyAl = al + 10*randn(size(al,1));
blurNoisyAl = conv2(noisyAl,uniformFilter,'same');
showIm(noisyAl + i*blurNoisyAl);
drawnow

% The noise is attenuated, but the blurring is disturbing.  We'll
% see how to do better than this later.

% One can also take a weighted average using a smoothly decaying
% window -- here we make a Gaussian weighted window:
%
gaussFilter = mkGaussian([9 9],4);
gaussFilter = gaussFilter / sum2(gaussFilter);
clf; showIm(gaussFilter);

% The weights must sum to 1 -- think about why this should be
% when computing a local average:
sum(sum(gaussFilter))

% Here we use the Gaussian window to blur an image:
%
blurAlGauss = upConv(al,gaussFilter);
showIm(al + i*blurAlGauss)
drawnow

% Let's add some noise to Al, and see what the blurring does:
noisyAl = al + 10*randn(size(al,1));
showIm(noisyAl + i*upConv(noisyAl,mkGaussian([7 7],2.5)));
drawnow

% The blurring looks about the same with the gaussian and the
% uniform operators, but it's not.  These two operators have some
% significant differences.  Consider what they do to 2 sinusoids
% with different wavelengths:
%
s1 = mkSine([64 64],4);
s2 = mkSine([64 64],8);
u1 = upConv(s1,uniformFilter);
u2 = upConv(s2,uniformFilter);
g1 = upConv(s1,gaussFilter);
g2 = upConv(s2,gaussFilter);

% Let's cut out a piece of the output from each filter and
% compare it to the input.  Of course we expect the main effect
% of the filter to be to attenuate the higher frequency (shorter
% wavelength) input signal:
%
yvals=8:1:14;
xvals=10:54;
clf, subplot(3,2,1)
showIm(s1(yvals,xvals),[-1,1]);
title('s1: period 4');
subplot(3,2,3)
showIm(g1(yvals,xvals),[-1,1]);
title('g1: gaussian blur');
subplot(3,2,5)
showIm(u1(yvals,xvals),[-1,1]);
title('u1: uniform blur');
subplot(3,2,2)
showIm(s2(yvals,xvals),[-1,1]);
title('s2: period 8');
subplot(3,2,4)
showIm(g2(yvals,xvals),[-1,1]);
title('g2: gaussian blur');
subplot(3,2,6)
showIm(u2(yvals,xvals),[-1,1]);
title('u2: uniform blur');
drawnow

% The figure shows the original waveform and the result of
% applying the two fitlers.  The higher frequency waveforms are
% indeed attenuated, but note something else.  The locations of
% crests (peaks) and troughs in the uniformFilter output are
% shifted with respect to the input in one but not the other
% frequency.  The edges in the image were shifted, which is NOT
% desirable.

% We can also see this by plotting slices through the outputs.
%
clf
subplot(4,1,1); plot(s1(20,10:50));
subplot(4,1,2); plot(u1(20,10:50));
subplot(4,1,3); plot(s2(20,10:50));
subplot(4,1,4); plot(u2(20,10:50));
drawnow

% The phase reversal is evident in the top two plots.  There is
% no phase reversal between the bottom two.  We'll explain this
% later when we learn how to analyse filters.

% Of course there are other types of operators for taking
% averages of images in each local region.  Many of these are
% nonlinear.  For example, if instead of computing the mean of
% each neighbourhood we computed the median, then the effective
% filter is quite nonlinear.

% A second type of local operator is a filter for enhancing edges
% in images.  They are often called band-pass filters, and often
% they only enhance edges within a specific range of
% orientations.  Here are two filters tuned to horizontal and
% vertical edges:
%
hfilt = [-0.107 0.0 0.107; 
         -0.245 0.0 0.245; 
	 -0.107 0.0 0.107];
vfilt = [-0.107 -0.245 -0.107; 
          0.0    0.0    0.0  ; 
	  0.107  0.245  0.107];
disc = mkDisc(64);
hResponse=upConv(disc,hfilt);
vResponse=upConv(disc,vfilt);
clf
showIm(hResponse + i*vResponse);
drawnow

% A third use of filters, is to decompose images into useful
% signal components. For example, we might wish to split an image
% into a lowpass (blurry) component and a highpass (edge enhanced)
% component.  Let's define filters to do this:
%
filterLP = [1 4 6 4 1]' * [1 4 6 4 1] / 256;
filterHP = zeros(5,5);
filterHP(3,3)=1;
filterHP=filterHP-filterLP;
   
lpAl = upConv(al,filterLP);
hpAl = upConv(al,filterHP);
clf, subplot(1,2,1)
showIm(lpAl);
title('Lowpass')
subplot(1,2,2)
showIm(hpAl);
title('Highpass')
drawnow

% These filters were designed so that we can actually reconstruct
% the original image from their outputs in a simple way; we just
% sum them:
%
imStats(lpAl+hpAl, al);

% This works because the filters are linear.  Together they form
% an invertible linear transformation.

%-----------------------------------------------------------------------
% 6. STATISTICAL IMAGE OPERATIONS 

% It's often the case that we want to know something about the
% structure of the images we are processing.  One way is to
% describe the statistics of the image.  This is often critical
% in image restoration (where the noise properties must be
% understood) and in image coding (where the statistical
% description of the image can be used to design efficient
% codes).  It is also used extensively in image analysis.  One
% example is the descrimination of different textural regions
% within an image (e.g. to distinguish water from forest in
% satellite images).

% For example we might want to know the mean and variance of the
% pixel values in an image.  Let's look at einstein
%
mean2(al)
var2(al)

% The function imStats prints both these out for you, as well as
% the maximum and minimum values in the image.
%
imStats(al)

% These numbers are called first-order statistics.  They describe
% properties of the probability density function (PDF) of the
% ensemble of image pixel values, independent of where those
% pixels came from in the image.  Related to the first-order PDF
% is the histogram.  It represents the number of pixels at which
% each intensity value occurs:
%
clf
hist(al(:),64)
drawnow

% The notation al(:) converts the 2D array/matrix holding the
% image into a vector.  If we don't convert the matrix into a
% vector before doing calling the hist function, it computes a
% histogram of each row, and plots each in a separate colour
% which is confusing.  Try it: 
% 
hist(al,64)
drawnow

% The pixel intensity PDF (probability distribution function)
% is simply the histogram divided by the number of pixels
% in the image, so that it integrates to 1 (as all PDFs must):
%
npixels = size(al,1) * size(al,2);
pdf = hist(al(:),64) / npixels;
plot(pdf)
drawnow

% The mean and variance (computed above) provide only crude
% information about the pdf.

% The image histogram has many uses, some of which are described
% in Castleman, Chapter 6 and 7.  For example, if one wants to
% know how large an area of the image (i.e., how many pixels) has
% an intensity (or brightness) above a certain value, then you
% can integrate the appropriate region of the histogram.  For
% example, let's store a histogram of al, with 64 bins, storing
% numbers of pixels in the jth bin with intensities between
% (j-1)*4 and (j*4-1).
%
histAl=hist(al(:),[0:4:255]);

% If we wanted 256 bins, then the jth element of histAl would
% store the number of pixels that have an intensity of
% (j-1):
%
histAl=hist(al(:),[0:1:255]);

% Why j and j-1 here?  Remember that indicies in matlab start at
% 1, but intensities in al go from 0 to 255.

% So how many pixels have intensities greater than or equal to
% 100?
%
sum(histAl(101:1:256))

% There are other ways to do this.  One is to create a binary image 
% (1 bpp) that is 0 when intensity is less than 100, and 1
% otherwise:
%
thresImage = al >= 100;
sum2(thresImage)

% I should have something more here on histogram manipulation,
% but I haven't done it yet.  For now, see Castleman, Chapter 6

% Another interesting descriptor of a distribution is the entropy
% of the pdf, defined by the integral of -pdf(x)*log2(pdf(x)).
%
entropy2(al)

% If you don't know about entropy, it's essentially a lower bound
% on the number of bits per pixel we need to store for lossless
% compression (in this case, if we only use first-order stats to our
% advantage).  Check out Shannon's book called 'A Mathematical
% Theory of Communication".... it's great.

% Above when we clipped al, we removed a lot of information. One
% would therefore expect the entropy to go down. 
%
entropy2(clip(al,100,150))

% In images, the values at nearby pixels are rarely independent
% of one another.  The first-order statistics do not capture this
% spatial dependence.  It is for this reason that we have second
% and higher-order statistics. 
 
% The power spectrum is an example of a set of second order
% statistics.

ws=size(al,1);
window = hamming(ws) * hamming(ws)';
clf; showIm(al.*window); 
drawnow

powerAl = powSpect(al.*window);
% Plot the log of the power spectrum, since the range of
% values are significant.  Also, we add a little bit to crop
% the plot from below.
lpAl = log(powerAl + 0.1);  
showIm(lpAl); 
drawnow

% Along a few of the scanlines:
clf;
plot([-ws/2:ws/2-1], lpAl(129,:), 'r'); axis([-ws/2 ws/2-1 -5 20]); hold on;
plot([-ws/2:ws/2-1], lpAl(100,:), 'g');
plot([-ws/2:ws/2-1], lpAl(64,:), 'b');
hold off;
drawnow

% Another second order statistic is 
% the image co-occurrence matrix.  For vectors, this is
% simply the cross-correlation.
% There are several ways to compute the cross-correlation.  One
% is to compute the cross-correlation integral directly. Another
% is to use the Fourier transform that we learn about later:
%

tmp = ifft2(powerAl);
xcorrAl = fftshift(tmp .* conj(tmp));
showIm(xcorrAl)
drawnow

% ... along the middle scanline
clf; plot([-ws/2:ws/2-1], xcorrAl(128,:)); 
axis([-ws/2 ws/2-1 0 max(max(xcorrAl))]);
drawnow;

%-----------------------------------------------------------------------
% 7. GEOMETRIC IMAGE OPERATIONS 

% The fourth basic type of image operations are manipulations of
% the spatial structure of the image intensity function.  We
% might wish to translate the image, or rotate the image, scale
% the image etc.  In fact, you have seen examples of "morphing"
% operations that animaters use to manipulate shapes (usually
% faces).  But geometric manipulations have many more
% applications than that. E.g.:
% 1. They are used to invert image distortions owing to
%    inaccuracies in the digitization camera or lens.  Some nice
%    examples are given in Castleman in Chapter 8.
% 2. Geometric operations are also central in algorithms that
%    measure motion from image sequences.  The deformations are
%    used to warp one frame to stabilize it with respect to a
%    neighbouring frame in the sequence.  When stabilized, you
%    have accounted for the motion.
% 3. A third use is the creatation of panorama's or mosaics from
%    a set of overlapping images.  This enables one to take a
%    bunch of frames that capture just a portion of the scene at
%    high resolution, rather than one low resolution picture.
%    The individual images then have to be aligned and warped to
%    be stitched together in the right way.

% There are some facilities built into matlab and others in the
% ise toolbox to enable geometric operations.  E.g., we can
% circularly shift an image to the right ('circular' because
% pixels that go off the right are redrawn on the left):
%
showIm(circularShift(al,40,0));
drawnow

% Or down:
%
showIm(circularShift(al,0,40));
drawnow

% Or both:
%
showIm(circularShift(al,40,40));
drawnow

% We can also rotate the image (e.g., by 25 degress):
%
showIm(imrotate(al,25));
drawnow

% There is also a function to apply linear deformations: These
% include rotation, as above, and also shears:
%
M = [1 0.5 0; 0 1 0];
w=warpAffine2(al,M);
index=isnan(w);
w(index) = 0;
clf; showIm(w);
drawnow

% We can also expand or contract the image:
%
M = [1 0 0; 0 1 0];
u=warpAffine2(al, 0.7*M); index=isnan(u); u(index) = 0;
v=warpAffine2(al, 1.4*M - [1 ; 1] * [0 0 25]); index=isnan(v); v(index) = 0;
displayImage(u+i*v);
drawnow

% Or we can do any combination of these linear warps, using
% matrix multiplication to combine the linear transformations.

% To create useful and interest deformations of image requires
% two ingredients: 
%   1) a method for creating the mapping of pixel locations in
%         one image to pixel locations in the other image.
%   2) a method for interpolating the image (i.e. finding 
%         image intensity values between image pixels.
