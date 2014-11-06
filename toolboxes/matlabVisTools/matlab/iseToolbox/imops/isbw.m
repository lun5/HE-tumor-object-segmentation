function y = isbw(x)
%ISBW True for binary images.
%   ISBW(A) returns 1 if A is a binary image and 0 otherwise.
%   An binary image is a logical matrix that contains only the
%   values 0 or 1. 
%
%   See also ISIND, ISGRAY.

%   Clay M. Thompson 2-25-93
%   Copyright (c) 1993-1996 by The MathWorks, Inc.
%   $Revision: 5.3 $  $Date: 1996/09/13 16:22:13 $

if islogical(x)
    y = ~any(x(:)~=0 & x(:)~=1);
else
    y = 0;
end

y = logical(double(y));




