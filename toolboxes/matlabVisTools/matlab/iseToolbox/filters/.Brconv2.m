function result = rconv2(im,filt)
% RCONV2: Convolution, with boundaries handled via
% reflection about the edge pixels (calls upConv).
%
% result = rconv2(image,filter)
%
% DJH, 8/96

result = upConv(im,filt,[1,1],[0,0],'reflect1');
