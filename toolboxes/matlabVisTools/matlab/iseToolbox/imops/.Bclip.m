function result = clip(im, min, max)
% clip(im, min, max)
%	Clips image values between min and max.
%	Min and max default to 0 and 1 respectively if not specified.
%

if ~exist('min')
  min=0;
end
if ~exist('max')
  max=1;
end

result=im;

index = find(im < min);
result(index) = min * ones(size(index));

index = find(im > max);
result(index) = max * ones(size(index));
