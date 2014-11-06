function ppmWrite(im, filename, comment)
% function ppmWrite(im, filename, comment)
% Write out an image to the named file in unsigned chars with header.
% The image format used is ppm for color images. If the file name is
% pix.ppm, and the color planes are in im(:,1:3), type:
%		ppmWrite(im, 'pix.ppm')
% comment is an optional comment string to be written in the ppm header.
% It can contain newline characters.  Default ''.

if ~exist('comment')
  comment = '';
end

fid = fopen(filename,'w');
[rows, cols, colours] = size(im);
if ~(colours == 3)
     error('Not an RGB image.');
end

imShuffle = uint8(zeros(3, cols, rows));
if ((min(im(:)) >= 0.0) & (max(im) <= 1.0))
  imShuffle(1, :, :) = round(reshape(im(:, :, 1), rows, cols)' * 255);
  imShuffle(2, :, :) = round(reshape(im(:, :, 2), rows, cols)' * 255);
  imShuffle(3, :, :) = round(reshape(im(:, :, 3), rows, cols)' * 255);
else
  imShuffle(1, :, :) = round(reshape(im(:, :, 1), rows, cols)' );
  imShuffle(2, :, :) = round(reshape(im(:, :, 2), rows, cols)' );
  imShuffle(3, :, :) = round(reshape(im(:, :, 3), rows, cols)' );
end

fprintf(fid,'P6\n');

%% Output comments
if length(comment)>0
  EOL = 10;  % End of line
  if comment(end) ~= EOL
    comment = [comment EOL];  % Add an end of line character to the comment.
  end
  lineEnds = findstr(EOL, comment);
  for j=1:length(lineEnds)
    if j == 1
      line = ['# ' comment(1:lineEnds(j))];
    else
      line = ['# ' comment((lineEnds(j-1)+1):lineEnds(j))];
    end
    fprintf(fid, line);
  end
end
%% Output dimension
fprintf(fid,'%d %d\n255\n', cols, rows);

fwrite(fid,imShuffle,'uint8');
fclose(fid);
