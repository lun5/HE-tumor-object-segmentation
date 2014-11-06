% RES = upConv(IM, FILT, EDGES, STEP, START, STOP, RES)
%
% Upsample matrix IM, followed by convolution with matrix FILT.  These
% arguments should be 1D or 2D matrices, and IM must be larger (in
% both dimensions) than FILT.  The origin of filt
% is assumed to be floor(size(filt)/2)+1.
%
% EDGES is a string determining boundary handling:
%    'circular' - Circular convolution
%    'reflect1' - Reflect about the edge pixels
%    'reflect2' - Reflect, doubling the edge pixels
%    'repeat'   - Repeat the edge pixels
%    'zero'     - Assume values of zero outside image boundary
%    'extend'   - Reflect and invert
%    'dont-compute' - Zero output when filter overhangs OUTPUT boundaries
%
% Upsampling factors are determined by STEP (optional, default=[1 1]),
% a 2-vector [y,x].
% 
% The window over which the convolution occurs is specfied by START 
% (optional, default=[1,1], and STOP (optional, default = 
% step .* (size(IM) + floor((start-1)./step))).
%
% RES is an optional result matrix.  The convolution result, once computed,
% will be destructively added into this matrix.  If this argument is passed,
% the result matrix will not be returned. DO NOT USE THIS ARGUMENT IF 
% YOU DO NOT UNDERSTAND WHAT THIS MEANS!!
% 
% NOTE: this operation corresponds to multiplication of a signal
% vector by a matrix whose columns contain copies of the time-reversed
% (or space-reversed) FILT shifted by multiples of STEP.  See corrDn.m
% for the operation corresponding to the transpose of this matrix.

% Eero Simoncelli, 6/96.  revised 2/97.
% Bug fix: ADJ, 10/01

function result = upConv(im,filt,edges, step,start,stop,res)

%% THIS CODE SHOULD NOT ACTUALLY BE USED! (MEX FILE IS CALLED INSTEAD)
%% The following code is just in case there is no mex file.  It is NOT
%% optimized.

fprintf(1,'Warning: You should compile the MEX code for "upConv", found in the MEX subdirectory.\n  It is much faster.\n');

%------------------------------------------------------------
%% OPTIONAL ARGS:

if (exist('edges') == 1) 
  if (strcmp(edges,'reflect1') ~= 1)
    warning('Using REFLECT1 edge-handling (use MEX code for other options).');
  end
end

if (exist('step') ~= 1)
  step = [1,1];
end	

if (exist('start') ~= 1)
  start = [1,1];
end	

% A multiple of step
if (exist('stop') ~= 1)
  stop = step .* (size(im) + floor((start-ones(size(start)))./step));
end	

if ( ceil((stop(1)+1-start(1)) / step(1)) ~= size(im,1) )
  error('Bad Y result dimension');
end
if ( ceil((stop(2)+1-start(2)) / step(2)) ~= size(im,2) )
  error('Bad X result dimension');
end

if (exist('res') ~= 1)
  res = zeros(size(im) .* step);
end	

%------------------------------------------------------------

tmp = zeros(size(res));
tmp(start(1):step(1):stop(1),start(2):step(2):stop(2)) = im;

result = zeros(size(res));
tmpSz = size(tmp);
fmid = ceil(size(filt)/2);
for i = 1:tmpSz(2)
  for j=1:tmpSz(1)
    for fi = 1:size(filt,2)
      ix = i-(fi-fmid(2));
      if ix < 1
        ix = 2-ix;
      end
      if ix >tmpSz(2)
        ix = 2*tmpSz(2) - ix;
      end  
      for fj = 1:size(filt,1)
        jy = j - (fj - fmid(1));
        if jy < 1
          jy = 2-jy;
        end
        if jy >tmpSz(1)
          jy = 2*tmpSz(1) - jy;
        end
        result(j, i) = result(j,i) + tmp(jy, ix) * filt(fj, fi);  
      end
    end
  end
end

res = result + res;
return;
 
