% IM = pgmread( FILENAME )
%
% Load a pgm image into a MatLab matrix.  
%   This format is accessible from the XV image browsing utility.
%   Only works for 8bit gray images (raw or ascii)

% Hany Farid, Spring '96.  Modified by Eero Simoncelli, 6/96.
% Modified by Allan Jepson, 9/99.
function im = pgmread( fname );

global SIZE_LIMIT;
if SIZE_LIMIT
 fprintf(1, 'SIZE_LIMIT recognized by pgmread\n');
end

[fid,msg] = fopen( fname, 'r' );
if (fid == -1)
  error(msg);
end

[pars type]= pnmReadHeader(fid);
if (pars==-1)
  fclose(fid);
  error([fname ': cannot parse pgm header']);
end

if ~(type == 'P2' | type == 'P5')
  fclose(fid);
  error([fname ': Not of type P2 or P5.']);
end

xdim = pars(1);
ydim = pars(2);
%%% Maximum pixel value
maxval = pars(3);
%fprintf(1, 'original image size: cols %d rows %d\n', xdim, ydim)
sz = xdim * ydim;

if SIZE_LIMIT
 if sz >= 16384
   ydim = floor(16384/xdim);
   sz = xdim * ydim;
   fprintf(1, 'truncated image size: cols %d rows %d\n', xdim, ydim)
 end
end

if (type == 'P2')
  [im,count]  = fscanf(fid,'%d',sz);
elseif (type == 'P5')

  count = 0;
  im = [];
  stat = fseek(fid, -sz, 'eof');
  if ~stat
    [im,count]  = fread(fid,sz,'uchar');
  end
  
  if (count ~= sz)
    fprintf(1,'Warning: File ended early! %s\n', fname);
    fprintf(1,'...Padding with %d zeros.\n', sz-count);
    im = [im ; zeros(sz-count,1)];
  end

else  %% This should not occur, we have already checked the type.
  fclose(fid);
  error([fname ': Not of type P2 or P5.']);
end

im = reshape( im, xdim, ydim )';
fclose(fid);

