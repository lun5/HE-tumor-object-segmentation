function result = gradientMag(f,edges)
% GRADIENTMAG: Gradient magnitude of an image, using upConv, and 5-tap
%              derivative filters designed by Eero Simoncelli.  
% 
%      result=gradientMag(im,edges)
% 
%      edges - any valid edge handler for upConv (default is 'circular').
%
% DJH '96

if ~exist('edges')
  edges='circular';
end

fx = dx(f,edges);
fy = dy(f,edges);
result = sqrt(fx.^2 + fy.^2);
