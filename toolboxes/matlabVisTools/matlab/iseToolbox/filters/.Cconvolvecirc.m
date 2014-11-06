% CONVOLVECIRC	circular convolution
%
%	convolvecirc(im,filter) performs a 1D or 2D circular convolution
%	of image im with filter.
%
%	convolvecirc(im,filter,[stepx,stepy],[startx,starty]) performs
%	convolution with subsampling step sizes (stepx,stepy) and starting
%	point (startx,starty) in the image.
%
%	C code written by Eero Simoncelli for Lisp based OBVIUS
%	converted to Matlab Mex file by Geoff Boynton on 11/11/94	
%
%	see also EXPANDCIRC


