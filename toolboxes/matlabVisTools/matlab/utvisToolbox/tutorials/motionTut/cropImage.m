function  [imC] = cropImage(im, cropBox)
  % imC = im(cropBox(1,2):cropBox(2,2), cropBox(1,1):cropBox(2,1));
  % Crop the image to the crop box.  Return an image of the size
  % dictated by cropBox.  It may be padded with values 128.
  % Here cropBox = [xMin, yMin; xMax, yMax];
    

  % Compute the size of the cropped image.
  sizeIm = size(im);
  sizeImC = cropBox(2,:) - cropBox(1,:) + 1;
  sizeImC = [sizeImC(2) sizeImC(1)];
  imC = 128*ones(sizeImC);
  
  % Intersect the cropBox with the image.
  cropBoxTrim = max( cropBox, 1);
  cropBoxTrim(1,:) = min( [cropBoxTrim(1,:); sizeIm(2) sizeIm(1)], [], 1);
  cropBoxTrim(2,:) = min( [cropBoxTrim(2,:); sizeIm(2) sizeIm(1)], [], 1);
  sizeImCT = cropBoxTrim(2,:) - cropBoxTrim(1,:) + 1;
  sizeImCT = [sizeImCT(2) sizeImCT(1)];

  % Copy the image data to the trimmed crop box with the cropped image.
  corner = (cropBox(1,:)-cropBoxTrim(1,:));
  if all(corner>=0) & corner(1) <= sizeImC(2) & corner(2) <= sizeImC(1)
    imC((corner(2)+(1:sizeImCT(1))), (corner(1)+(1:sizeImCT(2)))) = ...
        im(cropBoxTrim(1,2):cropBoxTrim(2,2), ...
           cropBoxTrim(1,1):cropBoxTrim(2,1));
  end

