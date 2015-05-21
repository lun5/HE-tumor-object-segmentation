function result = blur(im,level,filt,edges)
% BLUR: Blurs an image by blurring and subsampling repeatedly, followed by
% upsampling and blurring.
%
%      result=blur(im,level,filt,edges)
%
%      im - input image.
%      level - number of times to blur and subsample (default is 1).
%      filter - blurring 1d filter to be applied separably to the rows and
%               cols of im (default is (root2/16)*[1 4 6 4 1]).  Filter taps
%               must sum to root2.
%      edges - must be a valid edge handler for corrDn and upConv (default
%              is 'circular').
%
% DJH '96


if ~exist('level')
  level=1;
end

if ~exist('filt')
  filt = sqrt(2)*[0.0625 0.25 0.375 0.25 0.0625];
end
filt=filt(:)';

if ~exist('edges')
  edges='circular';
end

if level==0
  result=im;
else
  tmp1 = corrDn(im,filt,[1,2],[0,0],edges);
  tmp2 = corrDn(tmp1,filt',[2,1],[0,0],edges);
  tmp3 = blur(tmp2,level-1);
  tmp4 = upConv(tmp3,filt',[2,1],[0,0],edges,size(tmp1));
  result = upConv(tmp4,filt,[1,2],[0,0],edges,size(im));
end
