function [H, Sa] = linEstH3D(left, right, NUM_RESCALE)
  % [H, Sa] = linEstH3D(left, right, NUM_RESCALE)
  % Estimate the 3D homography matrix H from two 4 x n matrices of
  % corresponding points left and right.  
  % Here alpha(k) * left(:,k)  - (H * right(:,k)) apprx= 0 
  % NUM_RESCALE (default TRUE) uses Hartley's rescaling. Always use
  % rescaling, unless you wish to show how badly the un-normalized
  % algorithm works.
  % Returns H along with the singular values Sa of the 2nPts x 9 homogeneous
  % linear system for H.
    
  if nargin < 3
    NUM_RESCALE = 1;
  end
  
  nPts = size(left,2);
  if nPts < 5 | nPts ~= size(right,2)
    fprintf(2,'lineEstH: Innappropriate number of left and right points.');
    H = [];
    return;
  end
  
  if size(left,1) == 3
    left = [left; ones(1, nPts)];
  else % Normalize to pixel coords
    left = left./repmat(left(4,:), 4,1);
  end
  if size(right,1) == 3
    right = [right; ones(1, nPts)];
  else % Normalize to pixel coords
    right = right./repmat(right(4,:), 4,1);
  end
  
  imPts = cat(3, left, right);
  
  %% Rescale image data for numerical stability.
  if NUM_RESCALE
    Knum = repmat(eye(4), [1,1,2]);
    %%% Rescale for numerical stability
    mn = sum(imPts(1:3,:,:),2)/nPts;
    mns = reshape(mn, [3 1 2]);
    var = sum(sum((imPts(1:3,:,:)-repmat(mns, [1 nPts 1])).^2,2)/nPts, 1);
    %% Scale image points so that sum of variances of x and y = 2.
    scl = sqrt(3./var(:));
    %% Sanity: varScl =  var .* reshape(scl.^2, [1 1 2]) % Should be 3
    
    %% Scale so x and y variance is roughly 1, translate so image mean (x,y) is zero.
    Knum(1:3,4,:) = -mn;
    Knum(1:3,:,:) = Knum(1:3,:,:).*repmat(reshape(scl, [1 1 2]), [3, 4,1]);
    for kIm = 1:2
      imPts(:,:,kIm) = reshape(Knum(:,:,kIm),4,4) * imPts(:,:,kIm);
    end
    %% Sanity check
    % sum(imPts(1:3,:,:),2)/nPts  % Should be [0 0]'
    % sum(sum(imPts(1:3,:,:).^2,2)/nPts,1) % Should be 3.
  end

  %% Make constraint matrix A. 
  %% The matrix H satisfies: A h = 0, where h = (H_1,1; H_1,2; ... H_4,4).
  left = reshape(imPts(:,:,1), [4 nPts]);
  right = reshape(imPts(:,:,2), [4 nPts]);
  A = [];
  Id = eye(4);
  for k = 1:nPts
    [mx ix] = max(abs(imPts(:,k,1)));
    eix = Id(:,ix);
    C = kron(-left(ix,k)* eye(4), right(:,k)');
    C = C + kron(left(:,k)*eix', right(:,k)');
    %% Delete row ix
    Q = [];
    if ix > 1
      Q = C(1:(ix-1),:);
    end
    if ix<4
      Q = [Q; C((ix+1):4,:)];
    end
    A = [A; Q];
  end
  
  %% Factor A
  [Ua Sa Va] = svd(A); Sa = diag(Sa);
  
  %% Set H to be the right null vector of A, reshaped to a 3x3 matrix.
  H = reshape(Va(:,end), 4,4)';

  %% Undo the renormalization
  if NUM_RESCALE
    H = inv(reshape(Knum(:,:,1),4,4)) * H * reshape(Knum(:,:,2),4,4);
  end

  %% Modify H to make it norm 1.
  H = H / norm(H(:));


