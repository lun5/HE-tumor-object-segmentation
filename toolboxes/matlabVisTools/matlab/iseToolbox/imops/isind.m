function y = isind(x)
%ISIND True for indexed images.
%   ISIND(A) returns 1 if A is a nonlogical 2-D uint8 array or if
%   A is a 2-D double array containing integers that are greater
%   than or equal to 1.
%
%   See ISGRAY, ISBW.

%   Clay M. Thompson 2-25-93
%   revised by Chris Griffin 6-96
%   Copyright (c) 1993-1996 by The MathWorks, Inc.
%   $Revision: 5.4 $  $Date: 1996/09/18 21:57:32 $


if isa(x, 'uint8')
   if islogical(x)      % It's a binary image
      y = 0;
   elseif ndims(x)==2,  % Could be index or grayscale, we can't tell.
      y = 1;
   else                 % Most likely it's a mxnx3 RGB image
      y = 0;  
   end
else
    y = (min(x(:))>=1 & max(x(:))<=256) & all(abs(x(:)-floor(x(:)))<eps);
end    

y = logical(y);
