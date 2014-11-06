function [cropBox] = setCropBox(bbox, sigmaFilt, filtLen, sizeIm)
% [cropBox] = setCropBox(bbox, sigmaFilt, filtLen, sizeIm)
% Format of output crop box:
%    cropBox = [xBox(1) yBox(1); xBox(2) yBox(2)];
% Input:
%   bbox = [xBox(1) yBox(1); xBox(2) yBox(2)];
%   Enlarge bounding box bbox to have extra borders of width sigmaCrop.
%   Also, ensure the bounding box is at least the length of the filter
%%  (filtLen) wide and high,

  
%% Pad bounding box by sigmaFilt all around
cropBox = [(bbox(1,:)-1.0*sigmaFilt) ; (bbox(2,:)+1.0*sigmaFilt)];
cropBox = round(cropBox);

% Intersect the cropBox with the image.
cropBox = max( cropBox, 1);
cropBox(1,:) = min( [cropBox(1,:); sizeIm(2) sizeIm(1)], [], 1);
cropBox(2,:) = min( [cropBox(2,:); sizeIm(2) sizeIm(1)], [], 1);

% Ensure cropBox is at least filtLen in x direction
extra = ceil(filtLen - (cropBox(2,1) - cropBox(1,1) + 1));
if extra > 0 
  % Figure out space on left and right (spL spR)
  spL = min(cropBox(1,1) - 1, extra);
  spR = min(sizeIm(2) - cropBox(2,1), extra);
  if spL+spR < extra
    %% Not enough room
    cropBox = [];
  end
  % Figure out extra to add on left and right (exL exR)
  % Try to make the extra amount on each side symmetric 
  exL = ceil(extra/2);
  exR = extra - exL;
  if exL > spL
    exR = exR + (exL-spL);
    exL = spL;
  end
  if exR > spR
    exL = exL + (exR-spR);
    exR = spR;
  end
  
  cropBox(1,1) = cropBox(1,1)-exL;
  cropBox(2,1) = cropBox(2,1) + exR;
end

% Ensure cropBox is at least filtLen in y direction
extra = ceil(filtLen - (cropBox(2,2) - cropBox(1,2) + 1));
if extra > 0 
  % Figure out space on top and bottom (spL spR)
  spL = min(cropBox(1,2) - 1, extra);
  spR = min(sizeIm(1) - cropBox(2,2), extra);
  if spL+spR < extra
    %% Not enough room
    cropBox = [];
  end
  % Figure out extra to add on on top and bottom (exL exR)
  % Try to make the extra amount on each side symmetric 
  exL = ceil(extra/2);
  exR = extra - exL;
  if exL > spL
    exR = exR + (exL-spL);
    exL = spL;
  end
  if exR > spR
    exL = exL + (exR-spR);
    exR = spR;
  end
  
  cropBox(1,2) = cropBox(1,2)-exL;
  cropBox(2,2) = cropBox(2,2) + exR;
end

return;