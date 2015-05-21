function [F, Sa, Sf] = linEstF(left, right, NUM_RESCALE)
  % [F, Sa, Sf] = linEstF(left, right)
  % Estimate F matrix from 3 x n matrices of corresponding points left
  % and right.  
  % NUM_RESCALE (default TRUE) uses Hartley's rescaling. Always use
  % rescaling, unless you wish to show how badly the un-normalized
  % algorithm works.
  % Returns F,  the singular values Sa of the 9x9 linear system, and
  % Sf, the singular values of the approximate F matrix (before grabbing
  % a rank 2 approximation.
    
  
  if nargin < 3
    NUM_RESCALE = 1;
  end
  
  nPts = size(left,2);
  if nPts < 8 | nPts ~= size(right,2)
    fprintf(2, 'lineEstF: Innappropriate number of left and right points.');
    F = [];
    return;
  end
  
  if size(left,1) == 2
    left = [left; ones(1, nPts)];
  else % Normalize to pixel coords
    left = left./repmat(left(3,:), 3,1);
  end
  if size(right,1) == 2
    right = [right; ones(1, nPts)];
  else % Normalize to pixel coords
    right = right./repmat(right(3,:), 3,1);
  end
  
  imPts = cat(3, left, right);
  
  %% Rescale image data for numerical stability.
  if NUM_RESCALE
    Knum = repmat(eye(3), [1,1,2]);
    %%% Rescale for numerical stability
    mn = sum(imPts(1:2,:,:),2)/nPts;
    mns = reshape(mn, [2 1 2]);
    var = sum(sum((imPts(1:2,:,:)-repmat(mns, [1 nPts 1])).^2,2)/nPts, 1);
    %% Scale image points so that sum of variances of x and y = 2.
    scl = sqrt(2./var(:));
    %% Sanity: varScl =  var .* reshape(scl.^2, [1 1 2]); % Should be 2
    
    %% Scale so x and y variance is roughly 1, translate so image mean (x,y) is zero.
    Knum(1:2,3,:) = -mn;
    Knum(1:2,:,:) = Knum(1:2,:,:).*repmat(reshape(scl, [1 1 2]), [2, 3,1]);
    for kIm = 1:2
      imPts(:,:,kIm) = reshape(Knum(:,:,kIm),3,3) * imPts(:,:,kIm);
    end
    %% Sanity check
    % sum(imPts(1:2,:,:),2)/nPts  % Should be [0 0]'
    % sum(sum(imPts(1:2,:,:).^2,2)/nPts,1) % Should be 2.
  end

  %% Make constraint matrix A. 
  %% The matrix F satisfies: A f = 0, where f = (F_1,1; F_1,2; ... F_3,3).
  left = reshape(imPts(:,:,1), [3 nPts]);
  right = reshape(imPts(:,:,2), [3 nPts]);
  A = [(repmat(left(1,:)',1,3).* right') (repmat(left(2,:)',1,3).* right') ...
       (right')];

  %% Factor A
  [Ua Sa Va] = svd(A); Sa = diag(Sa);

  %% Set F to be the right null vector of A, reshaped to a 3x3 matrix.
  F = reshape(Va(:,end), 3,3)';

  %% Modify F to make it rank 2.
  [Uf Sf Vf] = svd(F); Sf = diag(Sf);
  Sf0 = Sf;
  Sf0(end) = 0.0;
  F = Uf * diag(Sf0) * Vf';

  %% Undo the renormalization
  if NUM_RESCALE
    F = reshape(Knum(:,:,1),3,3)' * F * reshape(Knum(:,:,2),3,3);
  end

