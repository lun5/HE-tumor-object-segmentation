% EXPANDCIRC	image expansion followed by circular convolution
%
%	expandcirc(im,filter) performs a 1D or 2D circular convolution
%	of image im with filter.
%
%	expandcirc(im,filter,[stepx,stepy],[startx,starty],[s_x,s_y]) performs
%	convolution with expansion step sizes (stepx,stepy), starting
%	point (startx,starty) in the image, and output image size
%	(s_x,s_y).  
%
%	C code written by Eero Simoncelli for Lisp-based OBVIUS
%	converted to Matlab Mex file by Geoff Boynton on 11/11/94	
%
%	see also CONVOLVECIRC GAUSSPYRAMID, LAPLACEPYRAMID

