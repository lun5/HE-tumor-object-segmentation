function  endPnts = cropLineInBox(n, c, cropBox)
%  endPnts = cropLineInBox(n, c, cropBox)
% For a line defined by n(1) x(1) + n(2) x(2) + c = 0, b
% and a box determined by cropBox = [x(1) y(1) x(2) y(2)] with
% x(1) <= x(2), y(1) <= y(2), return
% endPts = [x1 y1 x2 y2] where the line from x1,y1 to x2,y2
% is the segment of the original line cropped to be within
% crop box.  If there is no such line, then x1 is nan.

if (size(n,1) == 1) n = n'; end
%
% Check for sign change between successive vertices.
pts = [cropBox(1:2); cropBox(3:-1:2); cropBox(3:4); cropBox(1:3:4)];
err = pts * n + c;
cntPnts = 1;
endPnts = zeros(2,2);
for j= 1:4
  k = mod(j,4) + 1;
  if err(j) * err(k) < 0 
    % Sign change, find endpoint
    r = err(k)/(err(k) - err(j));
    endPnts(cntPnts, :) = pts(j,:) * r + (1-r) * pts(k,:);
    cntPnts = cntPnts+1;
  elseif err(k) == 0
    endPnts(cntPnts, :) = pts(k,:);
    cntPnts = cntPnts+1;
  end
  if cntPnts == 3
    break;
  end
end

if (cntPnts ~= 3)
  endPnts = endPnts * NaN;
end

return;

