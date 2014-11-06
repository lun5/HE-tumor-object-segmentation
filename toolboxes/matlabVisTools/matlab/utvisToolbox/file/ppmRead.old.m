function im = ppmRead(fname)
% function im = ppmRead(filename)
% Read an image from the named file in unsigned chars with header.
% The image format used is ppm for color images. If the file name
% is pix.ppm, type:
%		im = ppmRead('pix.ppm')
% Assumes a big endian byte order.
%

[fid,msg] = fopen( fname, 'r' );
if (fid == -1)
  error(msg);
end

[pars type]= pnmReadHeader(fid);
if (pars==-1)
  fclose(fid);
  error([fname ': cannot parse pgm header']);
end

if ~(type == 'P6')
  fclose(fid);
  error([fname ': Not of type P6.']);
end

xdim = pars(1);
ydim = pars(2);
%%% Maximum pixel value
maxval = pars(3);
fprintf(1, 'original image size: cols %d rows %d\n', xdim, ydim)
sz = xdim * ydim * 3;

[imShuffle,count]  = fread(fid,sz,'uint8');
fclose(fid);
if (count == sz)
  imShuffle = reshape( imShuffle, 3, xdim, ydim );
else
  fprintf(1,'Warning: File ended early!');
  imShuffle = reshape( [imShuffle ; zeros(sz-count,1)], 3, xdim, ydim);
end

im = zeros(ydim, xdim, 3);
im(:,:,1) = reshape(imShuffle(1, :, :), xdim, ydim)';
im(:,:,2) = reshape(imShuffle(2, :, :), xdim, ydim)';
im(:,:,3) = reshape(imShuffle(3, :, :), xdim, ydim)';


