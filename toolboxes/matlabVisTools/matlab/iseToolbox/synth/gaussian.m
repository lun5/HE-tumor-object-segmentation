function result = gaussian(y,x,sigma)
% GAUSSIAN: Example for how to use mkSyntheticImage
%
% DJH '96

rsqrd = (x.^2 + y.^2);
result = exp(-(rsqrd/sigma^2));
