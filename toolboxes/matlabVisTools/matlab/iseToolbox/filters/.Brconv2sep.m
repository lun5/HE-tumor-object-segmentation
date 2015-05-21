function result = rconv2sep(im,rowfilt,colfilt)
% RCONV2SEP: Separable, convolution using reflecting edge handler.
% 
%      result=rconv2sep(im,rowfilt,colfilt)
%
%      im - input image.
%      rowfilt - 1d filter applied to the rows
%      colfilt - 1d filter applied to the cols
%
% Example: foo=rconv2sep(im,[1 4 6 4 1],[-1 0 1]);
%
% DJH '96

rowfilt=rowfilt(:)';
colfilt=colfilt(:);

tmp = upConv(im,rowfilt,[1,1],[0,0],'reflect1');
result = upConv(tmp,colfilt,[1,1],[0,0],'reflect1');

