function [mb,nb] = bestblk(siz,k)
%BESTBLK Determine block size for block processing.
%   SIZ = BESTBLK([M N],K) returns, for an M-by-N image, the
%   optimal block size for block processing. K is a scalar
%   specifying the maximum row and column dimensions for the
%   block. SIZ is a 1-by-2 vector containing row and column
%   dimensions for the block.
%
%   [MB,NB] = BESTBLK([M N],K) returns the row and column
%   dimensions in MB and NB, respectively.
%
%   [...] = BESTBLK([M N]) uses the default value of 100 for K.
%
%   See also BLKPROC.

%   Clay M. Thompson
%   Copyright (c) 1993-96 by The MathWorks, Inc.
%   $Revision: 5.3 $  $Date: 1996/09/11 18:29:11 $

if nargin==1, k = 100; end % Default block size

%
% Find possible factors of siz that make good blocks
%

% Define acceptable block sizes
m = floor(k):-1:floor(min(ceil(siz(1)/10),ceil(k/2)));
n = floor(k):-1:floor(min(ceil(siz(2)/10),ceil(k/2)));

% Choose that largest acceptable block that has the minimum padding.
[dum,ndx] = min(ceil(siz(1)./m).*m-siz(1)); blk(1) = m(ndx);
[dum,ndx] = min(ceil(siz(2)./n).*n-siz(2)); blk(2) = n(ndx);

if nargout==2,
  mb = blk(1); nb = blk(2);
else
  mb = blk;
end
