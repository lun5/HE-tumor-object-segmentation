%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SamplingTutorial
%%%   David Heeger, Eero Simoncelli, and Patrick Teo 6/96.  
%%%   Based on OBVIUS tutorial by David Heeger and Eero Simoncelli.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure there are no name conflicts with previous work
% (First save any previous results you wish to keep.)
clear;

% Subsampling is best thought of in two steps.  The first step is
% multiplication by an impulse train.  The second step is to get
% rid of the zero values in between impulses.

% Start by defining a range for a discrete signal
%
x_range = [0:1:63];

% Make an impulse train:
%
imp_train = zeros(1,64);
select = find(rem(x_range,4)==0);
imp_train(select) = ones(size(select));

clf, subplot(2, 1, 1)
bar(x_range, imp_train, 'r');
axis([0 63 0 1.2]);
title('Impluse train')

% The Fourier transform of an impulse train is another impulse
% train.  In most cases, the spacing of the impulse trains will
% be different.
%
f_range = [-32:1:31]
mag_ft_train = abs(fftshift(fft(fftshift(imp_train))));
subplot(2, 1, 2)
bar(f_range, mag_ft_train, 'g');
axis([-32 31 0 18]);
title('Magnitude of FFT of impluse train')
drawnow

% Multiplication in space yields convolution in frequency.  Thus,
% multiplying a signal by an impulse train yields a frequency
% domain representation with periodic replicas.
%
clf
gauss_12 = exp(-(((x_range-32).^2)/(12^2)));
subplot(2,2,1);
plot([-32:31], gauss_12, 'r-');
axis([-32 31 0 1]);  title('continuous function');

% Sample the gaussian by multiplying with the impulse train:
%
samp_gauss_12 = imp_train .* gauss_12;
subplot(2,2,3);
bar(x_range, samp_gauss_12, 'r');
axis([0 63 0 1]); title('sampled function');

% Compute the Fourier transform of the original gaussian:
%
mag_ft_gauss_12 = abs(fftshift(fft(fftshift(gauss_12))));
subplot(2,2,2);
plot([-32:31], mag_ft_gauss_12, 'g-');
axis([-32 31 0 24]); title('DFT of continuous function')

% Compute the Fourier transform of the sampled gaussian.  Notice
% the replicas.
%
mag_ft_samp_gauss_12 = abs(fftshift(fft(fftshift(samp_gauss_12))));
subplot(2,2,4);
plot([-32:31], mag_ft_samp_gauss_12, 'g-');
axis([-32 31 0 6]); title('DFT of sampled function')
drawnow

% If we make the Gaussian narrower in the space domain, then it
% becomes broader in the frequency domain.  If we make it narrow
% enough in space, then the frequency domain replicas start to
% overlap.

% Narrow gaussian:
%
clf
gauss_4 = exp(-(((x_range-32).^2)/(4^2)));
subplot(2,2,1);
plot([-32:31], gauss_4, 'r-');
axis([-32 31 0 1]);title('continuous function');

% Sample the narrow gaussian:
%
samp_gauss_4 = imp_train .* gauss_4;
subplot(2,2,3);
bar(x_range, samp_gauss_4, 'r');
axis([0 63 0 1]);title('sampled function');

% Fourier transform of narrow gaussian.  Notice the result is a
% wider gaussian.
%
mag_ft_gauss_4 = abs(fftshift(fft(fftshift(gauss_4))));
subplot(2,2,2);
plot([-32:31], mag_ft_gauss_4, 'g-');
axis([-32 31 0 8]); title('DFT of cts function')


% The Fourier transform of the sampled gaussian is so wide that
% it overlaps.
%
mag_ft_samp_gauss_4 = abs(fftshift(fft(fftshift(samp_gauss_4))));
subplot(2,2,4);
plot([-32:31], mag_ft_samp_gauss_4, 'g-');
axis([-32 31 0 2]);title('DFT of sampled function')
drawnow

% In this example, the replicas overlap each other.  This situation is
% called "aliasing". It is called aliasing because the sampled signal
% is masquerading as a different signal.  For example, let's look at
% two sinusoids that alias to one another when sampled.

% Low-frequency sinusoid:
%
sine_4 = -sin(2*pi*(4/64)*x_range);
clf, subplot(2,2,1);
plot([-32:31], sine_4, 'r-');
axis([-32 31 -1.2 1.2]);
title('continuous low freq')

% High-frequency sinusoid:
%
sine_12 = sin(2*pi*(12/64)*x_range);
subplot(2,2,2);
plot([-32:31], sine_12, 'y-');
axis([-32 31 -1.2 1.2]);
title('continuous high freq')

% Sample the low-frequency sinusoid:
%
samp_sine_4 = imp_train .* sine_4;
subplot(2,2,3);
bar(x_range, samp_sine_4, 'r');
axis([0 63 -1.2 1.2]);
title('sampled low freq')

% Sample the high-frequency sinusoid:
%
samp_sine_12 = imp_train .* sine_12;
subplot(2,2,4);
bar(x_range, samp_sine_12, 'y');
axis([0 63 -1.2 1.2]);
title('sampled high freq')
drawnow

% Even though sine-4 and sine-12 are distinct signals, their
% sampled versions are indistinguishable.  Sampled-sine-4 is an
% "alias" for sampled-sine-12.

% Nyquist Theorem states that to avoid aliasing a bandlimited
% signal should be sampled at twice the highest frequency.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction:

% If there is no aliasing, we can reconstruct the original signal
% from the sampled signal.  One easy way to do this is in the
% frequency domain, by getting rid of the replicas, and then
% using the inverse Fourier transform.  Why is that the right
% thing to do?  Well, let's go back and compare the Fourier
% transform of the original signal to that of the sampled signal.

clf;

% Fourier transform of gaussian:
%
subplot(2, 1, 1), plot([-32:31], mag_ft_gauss_12, 'r-');
axis([-32 31 0 24]);
title('Fourier transform of gaussian')

% Fourier transform of sampled gaussian:
%
subplot(2, 1, 2)
plot([-32:31], mag_ft_samp_gauss_12, 'g-');
axis([-32 31 0 6]);
title('Fourier transform of sampled gaussian')

% The main difference between the two are that there are replicas
% in the second. There's also a scale factor of 4 (note that we
% subsampled by a factor of 4), but that'll be easy to fix.  To
% get rid of the replicas in the frequency domain, we multiply by
% a box (also known as an "ideal" low-pass filter).
%
box = zeros(1,64);
select = find(x_range>24 & x_range<40);
box(select) = ones(size(select));
subplot(2, 1, 2), hold on
plot(x_range-32, box, 'c-');
hold off
drawnow

% Use the box to reconstruct to the (non-aliased) Gaussian: (by
% mulitplying the two together in the frequency domain. Note the
% scale factor of 4 to account for the subsampling by a factor of
% 4.)
%
clf
subplot(2,1,1);
plot([-32:31], gauss_12, 'r-');
axis([-32 31 0 1]);  title('original function');

recon_ft_gauss_12 = box .* fftshift(fft(fftshift(samp_gauss_12)));
recon_gauss_12 = 4.*real(fftshift(ifft(fftshift(recon_ft_gauss_12))));
subplot(2,1,2);
plot([-32:31], recon_gauss_12, 'r-');
axis([-32 31 0 1]); title('reconstruction');
drawnow

% Should be zero (the mean squared error between original and reconstructed
% version):
%
mse(recon_gauss_12, gauss_12)

% However, for the aliased Gaussian it doesn't work.  The box
% filter pulls out one section of the frequency domain.  In this
% case, the section that is pulled out is already contaminated by
% parts of the neighboring replicas.
%
clf
subplot(2,1,1);
plot([-32:31], gauss_4, 'r-');
axis([-32 31 -0.2 1.2]);  title('original function');

recon_ft_gauss_4 = box .* fftshift(fft(fftshift(samp_gauss_4)));
recon_gauss_4 = 4.*real(fftshift(ifft(fftshift(recon_ft_gauss_4))));
subplot(2,1,2);
plot([-32:31], recon_gauss_4, 'r-');
axis([-32 31 -0.2 1.2]); title('reconstruction');
drawnow

% Multiplying by the box in the frequency domain is the same as
% convolving with its inverse Fourier transform in the space
% domain.  The inverse Fourier transform of a box is called a
% sinc filter:
%
clf
sinc = real(fftshift(ifft(fftshift(box))));
plot([-32:31], sinc, 'r-');
axis([-32 31 -0.1 0.3]); title('sinc function')
drawnow

% Now you see why it's easier to do reconstruction in the
% frequency domain.  Convolution with such a big filter would
% take a long time.

% We don't always need to use sinc interpolation to reconstruct
% the original signal.  If there is room to spare between the
% replicas then we can use a filter with a more gradual fall-off.
% For example, let's reconstruct the sine-4 from sampled-sine-4.
% First let's use sinc interpolation (as before):
%
recon_ft_sine_4 = box .* fftshift(fft(fftshift(samp_sine_4)));
recon_sine_4 = 4.*real(fftshift(ifft(fftshift(recon_ft_sine_4))));
clf
subplot(2,1,1)
plot([-32:31], sine_4, 'r-');
axis([-32 31 -1.2 1.2]); title('original function')
subplot(2,1,2)
plot([-32:31], recon_sine_4, 'r-');
axis([-32 31 -1.2 1.2]); title('reconstruction')
drawnow

% Should be zero:
mse(recon_sine_4, sine_4)

% Now let's use a different filter, that falls off smoothly on
% both sides:
%
cos_box = zeros(1,64);
for x = 0:63
  if (x>20 & x<28)
    cos_box(x+1) = 0.5*(1+cos(pi*(x-28)/8));
  elseif (x>27 & x<37)
    cos_box(x+1) = 1;
  elseif (x>36 & x<44)
    cos_box(x+1) = 0.5*(1+cos(pi*(x-20)/8));
  end
end
clf
plot(x_range, cos_box, 'c-');
axis([0 63 -0.2 1.2]);
title('smooth box')
drawnow

% Reconstruct sampled sinusoid:
%
new_recon_ft_sine_4 = cos_box .* fftshift(fft(fftshift(samp_sine_4)));
new_recon_sine_4 = 4.*real(fftshift(ifft(fftshift(new_recon_ft_sine_4))));
clf
subplot(2,1,1)
plot([-32:31], sine_4, 'r-');
axis([-32 31 -1.2 1.2]); title('original function')
subplot(2,1,2)
plot([-32:31], new_recon_sine_4, 'r-');
axis([-32 31 -1.2 1.2]); title('reconstruction')

% Should be small:
mse(new_recon_sine_4, sine_4)

% Now let's look at the impulse response of this cos_box filter:
%
cos_filter = real(fftshift(ifft(fftshift(cos_box))));
clf
plot([-32:31], cos_filter, 'r-');
axis([-32 31 -0.1 0.3]);
title('Impulse response of smooth box')

% There are two points to be made here.  First, there is in
% general more than one way to reconstruct a subsampled signal.
% We have some freedom in choosing the interpolation filter.
% Second, the impulse response of the cos-box filter much more
% compact (in space) than the sinc filter, so in some cases we
% might decide to use convolution to reconstruct rather than
% reconstructing in the Fourier domain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% So far, we have only dealt with the first step of subsampling,
% multiplying by an impulse train.  The second step is to get rid
% of the zero values in between impulses.  Let's do this for the
% gaussian (the one that was not aliased).

% Subsampled gaussian with zeros removed:
%
sub_gauss_12 = samp_gauss_12([1:4:64]);
clf, 
subplot(2, 1, 1), plot([-8:7], sub_gauss_12, 'r-');
axis([-8 7 0 1]); title('subsampled gaussian')

% Fourier transform of subsampled gaussian:
%
mag_ft_sub_gauss_12 = ...
    abs(fftshift(fft(fftshift(sub_gauss_12))));
subplot(2, 1, 2), plot([-8:7], mag_ft_sub_gauss_12, 'g-');
axis([-8 7 0 6]); title('Magnitude of DFT');
drawnow

% Notice that this pulls out one cycle of the replicated
% frequency domain.

% It is easy to reconstruct the original signal from the fully
% subsampled signal, just by padding out the frequency domain with
% zeroes.
%
recon_ft_gauss_12 = zeros(1,64);
recon_ft_gauss_12([25:40]) = fftshift(fft(fftshift(sub_gauss_12)));
recon_gauss_12 = 4.*real(fftshift(ifft(fftshift(recon_ft_gauss_12))));
clf
subplot(2,1,1)
plot([-32:31], gauss_12, 'r-');
axis([-32 31 0 1]); title('original function')
subplot(2,1,2)
plot([-32:31], recon_gauss_12, 'r-');
axis([-32 31 0 1]); title('reconstruction')

% Note also the multiplication by 4.  The replicas of the spectrum are
% actually reduced in amplitude by a factor equal to the subsampling
% factor. (Recall the 1/N in the Matlab version of the inverse Fourier
% transform.)  In order to recover the signal, we have to multiply by this
% factor.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Why is aliasing such a bad thing?  After all, we are not
% necessarily interested in reconstructing the original signal.
% Here's an example that demonstrates how serious a problem
% aliasing can be for motion analysis.  This example also
% demonstrates how to reduce aliasing by using a pre-filter, that
% is, applying a low-pass filter before subsampling.

% Make a slightly tilted line:
%
nn_range = [-63.5:1:63.5];
[nn_rangeX,nn_rangeY] = meshgrid(nn_range);
theta = pi/16;
ramp = cos(theta).*nn_rangeX-sin(theta).*nn_rangeY;
select = find(abs(ramp)<1.5);
myLine = zeros(128,128);
myLine(select) = ones(size(select));

% Make an image sequence of the line moving slowly to the right:
%
line_seq = zeros(128*128,10);
for frame = 1:10
  fprintf(1, 'Creating image %d:\n', frame);
  shifted_line = circularShift(myLine, 2*(frame-1), 0);
  line_seq(:, frame) = shifted_line(:);
end

% Make a matlab movie and display it:
%
clf
showIm(myLine,'auto1','auto',0);
line_movie=moviein(10);
for frame=1:10
  showIm(reshape(line_seq(:,frame),128,128),'auto1','auto',0);
  line_movie(:,frame)=getframe;
end
movie(line_movie,5);

% Now we subsample the images in the sequence.  
%
sub_line_seq=zeros(32*32,10);
for frame=1:10
  tmp=reshape(line_seq(:,frame),128,128);
  tmp=tmp([1:4:128],[1:4:128]);
  sub_line_seq(:,frame) = tmp(:);
end

% Make movie and display it
%
showIm(reshape(sub_line_seq(:,1),32,32),'auto','auto',0);
sub_line_movie = moviein(10);
sub_line_movie=moviein(10);
for frame=1:10
  showIm(reshape(sub_line_seq(:,frame),32,32),'auto','auto',0);
  sub_line_movie(:,frame)=getframe;
end
movie(sub_line_movie,5);

% Because of the aliasing, the subsampled sequence appears to
% move UP!

% Let's try prefiltering the sequence before sampling.  When
% viewing this blurred sequence, notice that there is still a
% little bit of upward motion.  That is because the original
% line-seq images were already slightly aliased.

% Blur image sequence:
%
blur_line_seq=zeros(size(line_seq));
for frame = 1:10
  fprintf(1, 'Blurring frame %d:\n', frame);
  tmp=blur(reshape(line_seq(:,frame),128,128),2);
  blur_line_seq(:,frame)=tmp(:);
end

% Now we subsample the images in the sequence.  
%
sub_blur_line_seq=zeros(32*32,10);
for frame=1:10
  tmp=reshape(blur_line_seq(:,frame),128,128);
  tmp=tmp([1:4:128],[1:4:128]);
  sub_blur_line_seq(:,frame) = tmp(:);
end

% Make movie and display it:
%
showIm(reshape(sub_blur_line_seq(:,1),32,32),'auto','auto',0);
sub_blur_line_movie = moviein(10);
sub_blur_line_movie=moviein(10);
for frame=1:10
  showIm(reshape(sub_blur_line_seq(:,frame),32,32),'auto','auto',0);
  sub_blur_line_movie(:,frame)=getframe;
end
movie(sub_blur_line_movie,5);

% When we subsample the blurred sequence, there isn't much
% aliasing.  The line appears to move mostly rightward.

% What is going on?  It is the high frequency components that are
% aliasing.  The low-pass (blurring) filter attenuates the high
% frequency components, thereby reducing the impact of the
% aliasing.

% Not any filter will do.  Since we are subsampling by a factor
% of 4, you might think that averaging over a 4x4 patch will fix
% things.  NOT!  Try it.  

% Why doesn't a 4x4 local averaging filter work?  Look at its
% frequency response:
%
clf;
one_d_box = 0.25*ones(1,4);
impulse = mkImpulse([1,64]);
impulse_response = cconv2(impulse, one_d_box);
mag_freq_resp = abs(fftshift(fft(fftshift(impulse_response))));
plot([-32:31], mag_freq_resp, 'r-');
axis([-32 31 0 1]);
drawnow

% Since we are subsampling by a factor of 4, we want the
% pre-filter to pull out one-fourth of the frequencies (the low
% freqs), attenuating the others to zero (or nearly zero).
% Clearly, this filter doesn't do that.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reconstruct uniform samples from irregular samples:

% Given irregular samples of an appropriately bandlimited signal,
% you can usually compute what the uniform samples would have
% been (had the sampling been uniform).  The basic idea is that
% for an appropriately bandlimited signal, you can reconstruct
% the continuous signal from the irregular samples and then
% resample it uniformly.  Better yet, you can combine those two
% steps and go directly from the irregular samples to the uniform
% samples.  Here, we go through an example of how to do this.

% First, we construct some sampling matrices.  The first is a
% uniform subsampling that takes every 4th sample.  The second is
% a "jittered" subsampling in which each sample is from one of 4
% positions.

% Uniform sampling matrix:
%
reg_sampling_mat = zeros(16,64);
for s = 0:15
  reg_sampling_mat(s+1,s*4+1) = 1;
end
clf, subplot(2, 1, 1)
showIm(reg_sampling_mat); 
title('Matrix corresponding to uniform sampling')

% Irregular sampling matrix:
%
irreg_sampling_mat = zeros(16,64);
for s = 0:15
  irreg_sampling_mat(s+1,s*4+floor(rand([1,1])*4)+1) = 1;
end
subplot(2, 1, 2)
showIm(irreg_sampling_mat); 
title('Matrix corresponding to irregular sampling')
drawnow

% Next, we need to explain what we mean by an appropriately
% bandlimited signal.  For purposes of this example we use the
% Hartley basis set, that is a set of cas=sin+cos functions.  The
% matrix for Hartley transform is:
%
Hartley_mat = zeros(64,64);
nn_range = [0:63];
for k = 0:63
  Hartley_mat(k+1,:) = 1/sqrt(64).*(cos(2*pi/64*k.*nn_range)+...
      sin(2*pi/64*k.*nn_range));
end
clf, showIm(Hartley_mat); 
title('Matrix corresponding to Hartley transform basis functions')
drawnow

% The Hartley transform is an orthogonal transform (i.e., it is
% its own inverse):
%
showIm(Hartley_mat' * Hartley_mat); 
title('Should be identity matrix')

% For a bandlimited signal, we need only some of the basis
% functions to represent the signal.  Here are the first quarter
% of the basis functions:
%
Hartley_prime_mat = zeros(16,64);
nn_range = [0:63];
for k = 0:15
  Hartley_prime_mat(k+1,:) = 1/sqrt(64).*(cos(2*pi/64*k.*nn_range)+...
      sin(2*pi/64*k.*nn_range));
end
showIm(Hartley_prime_mat); 
title('Matrix corresponding to subset of Hartley basis functions')
drawnow

% Note that this modified Hartley transform is no longer
% generally invertible:
%
clf
showIm(Hartley_prime_mat'*Hartley_prime_mat);
title('Not identity matrix')

% Even so, there are some signals that can be adequately
% represented by the modified Hartley basis set.  We can
% construct such a signal in the transform domain:
%
H_transform = rand([16,1])*2-1;
cont_signal = Hartley_prime_mat' * H_transform;
clf
plot([0:63], cont_signal, 'r-');
axis([0 63 min(cont_signal) max(cont_signal)]);
title('Input signal: random combination of Hartley subset')

% Now check that this signal can be reconstructed from its
% Hartley transform.  Should be zero:
%
mse(Hartley_prime_mat'*(Hartley_prime_mat*cont_signal),...
    cont_signal)

% Finally, we construct a matrix that converts directly from the
% irregular samples to the regular samples.
%
irreg_to_reg_mat = (reg_sampling_mat*Hartley_prime_mat') * ...
    inv(irreg_sampling_mat*Hartley_prime_mat');

% There is a potential problem in the above calculation if the
% resulting matrix is not invertible.  This might happen
% depending exactly on where the irregularly spaced samples are.
% Check that the determinant is not too small:
%
sqrt(abs(det(irreg_sampling_mat*Hartley_prime_mat')))

% If this value is smaller than 1e-6 or so, then recompute 
% the irregular-sampling-matrix, and recompute the 
% irreg-to-reg-matrix.

% Why is this matrix singular for certain choices of
% irregular-sampling-matrix?  What is the condition on
% irregular-sampling-matrix and Hartley-prime-mat that would
% guarantee invertibility?

% Let's see how it works.  Sample the "continuous" signal:
%
reg_sampled = reg_sampling_mat * cont_signal;
clf, subplot(3, 1, 1), hold on
plot([0:15], reg_sampled, 'r-');
title('uniform samples')

irreg_sampled = irreg_sampling_mat * cont_signal;
subplot(3, 1, 2)
plot([0:15], irreg_sampled, 'b-');
title('irregular samples')
drawnow

% Notice how different those subsampled signals are.  Now,
% construct the regularly subsampled signal from the irregularly
% subsampled signal:
%
recon_reg_sampled = irreg_to_reg_mat * irreg_sampled;
subplot(3, 1, 3)
plot([0:15], recon_reg_sampled, 'b-'); 
title('reconstruction from irregular samples')

% ... superimposed on original (top) figure
subplot(3, 1, 1)
plot([0:15], recon_reg_sampled, 'bo'); 
hold off

% Should be zero:
mse(recon_reg_sampled, reg_sampled)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Local Variables:
%%% buffer-read-only: t 
%%% End:

