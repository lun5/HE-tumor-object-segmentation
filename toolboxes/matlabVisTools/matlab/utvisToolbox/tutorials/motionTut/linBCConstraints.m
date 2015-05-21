function [C] = linBCConstraints(im1, im0, gradIm0, idIm0, xMesh, yMesh, gTol2)
%% [C] = linBCConstraints(im1, im0, gradIm0, idIm0, gTol2)
%% Compute linearized brightness constancy contraints.
%%
%% Input: 
%%   Successive image frames, im0, im1.  
%%   gradIm0 is the gradient image of size [size(im0) 2].  
%%   idIm0 the interior of the bounding box in frame 0.
%%   xMesh and yMesh are the x and y coordinate images of size(im0).
%%   gTol2 minimum gradient amplitude to use for estimating motion
%%         constraints
  
  diffIm = im1 - im0;
  nrm2Im = sum(gradIm0.^2,3);
  idx = idIm0(:) & (nrm2Im(:) > gTol2);

  gradIm0 = reshape(gradIm0, prod(size(im1)), 2);
  nC = sum(idx);
  if nC > 0
    C = zeros(5, nC);
    C(1:2, :) = gradIm0(idx,1:2)';
    C(3, :) = diffIm(idx)';
    C(4:5, :) = [xMesh(idx) yMesh(idx)]';
  else
    C =[];
  end
  
  return;