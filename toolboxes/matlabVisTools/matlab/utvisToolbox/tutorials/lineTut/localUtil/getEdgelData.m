function [p, theta] = getEdgelData(edgelImage, dirImage, cropBox)
% [p, theta] = getEdgelData(edgelImage, dirImage, cropBox)
% Get the edgel data from within the cropBox, given
%    - edgelImage  - binary image marking edgel locations
%    - dirIm       - theta/pi for the edgel normal.
% Optional:
%    - cropBox     - [x0 y0; x1 y1] for top left (x,y) = (x0,y0)
%                    and bottom right corners of crop box.
%                  - [] or absent, use whole image.
% This returns:
%   - p      an M by 2 array for the x,y image positions
%   - theta  an M by 1 vector for the edgel orientations
%            The edgel normal is (sin(theta), cos(theta))
% ADJ Fall 00.

if (nargin < 3)
  cropBox = [];
end

if size(cropBox,1) > 0
  cropBox(cropBox < 1) = 1;
  if cropBox(1) > size(edgelImage,2)
     cropBox(1) = size(edgelImage,2);
  end
  if cropBox(3) > size(edgelImage,2)
     cropBox(3) = size(edgelImage,2);
  end
  if cropBox(2) > size(edgelImage,1)
     cropBox(2) = size(edgelImage,1);
  end
  if cropBox(4) > size(edgelImage,1)
     cropBox(4) = size(edgelImage,1);
  end
else %% Default cropBox
  cropBox = [1 1 size(edgelImage, 2) size(edgelImage,1)];
end

cropEdgels = edgelImage(cropBox(2):cropBox(4), cropBox(1):cropBox(3));
cropDir = dirImage(cropBox(2):cropBox(4), cropBox(1):cropBox(3));
[x, y] = meshgrid(cropBox(1):cropBox(3), cropBox(2):cropBox(4));
edgelIndices = cropEdgels > 0;
p = [x(edgelIndices) y(edgelIndices)];
theta = cropDir(edgelIndices) * pi;
