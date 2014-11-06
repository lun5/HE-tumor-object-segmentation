function [aff,outliers,robustError] = fit_robust_affine_transform( p0, p1, wght, est_outlier_rate_or_init_aff )

% [aff] = fit_robust_affine_transform( p0, p1 )
%
% Robustly compute the affine transformation that takes points
% p0 to points p1 using RANSAC and reweigted least squares.
%
% Input:
% p0 - initial points stores as an 2xn matrix.
% p1 - points after transformation.
% wght - weight of each constraint.
% est_outlier_rate_or_init_aff - the estimated outlier rate or
%   the initial affinie transform.  It an outlier rate is specified,
%   RANSAC is used to get an initial guess.  If an affine transform is
%   specified, it is used instead.
%
% Output:
% aff - affine tranformation embeded in a 3x3 matrix.
% outliers - indices of the outliers to the estimated
%   transformation.
% robustError - robust error of the fit.
%
% Thomas F. El-Maraghi
% May 2004

if ~exist('wght')
   wght = ones(1,size(p0,2));
end

if ~exist('est_outlier_rate_or_init_aff')
   est_outlier_rate = 0.5;
else
   if length(est_outlier_rate_or_init_aff(:)) == 1 
      est_outlier_rate = est_outlier_rate_or_init_aff;
   else
      init_aff = est_outlier_rate_or_init_aff;
   end
end

% Check that there are enough constraints
n = size(p0,2);
if n < 3
   error( 'Too few contraints to fit affine tranform.' );
end

if ~exist('init_aff')
   % Perform RANSAC to generate initial guess for the
   % affine transformation
   num_guesses = ceil(1.0/(1.0-est_outlier_rate)^3);
   best_seed = zeros(3,1);
   best_median = Inf;
   best_aff = zeros(3,1);
   for k = 1:num_guesses
      % Select 3 random constraints
      c = [ceil(rand * n)];
      for j = 2:3
         k = ceil(rand * n);
         while any(c == k)
            k = ceil(rand * n);
         end
         c = [c; k];
      end   
      
      % Build the least squares matrices   
      A = [];
      b = [];
      for j = 1:3
         x = p0(1,c(j));
         y = p0(2,c(j));
         A = [A; wght(j)*[x 0 y 0 1 0; 0 x 0 y 0 1]];
         b = [b; wght(j)*p1(:,c(j))];
      end
      
      % Solve for the affine tranform
      aff = pinv(A)*b;
      aff = [aff(1:2) aff(3:4) aff(5:6); 0 0 1];
      
      % Compute the residuals
      pts = aff * [p0; ones(1,n)];
      resid = sqrt(sum((pts(1:2,:) - p1).^2,1));
      
      % Save the transformation that has the lowest median residual
      if median(resid) < best_median
         best_seed = c;
         best_median = median(resid);
         best_aff = aff;
      end
   end
   aff = best_aff;
else
   aff = init_aff;
end
   
% Perform reweighted least squares to improve the estimate of
% the affine transform.
for k = 1:10
   % Compute the residuals
   pts = aff * [p0; ones(1,n)];
   resid = sqrt(sum((pts(1:2,:) - p1).^2,1))';
      
   % Compute the contraint weights
   sigma = 1.4826*median(resid);
   dev = resid/sigma;
   resid_wght = exp(-(dev.^2)/2)/(sqrt(2*pi)*sigma);
   resid_wght( find(dev>5.0) ) = 0;
   
   % Build the least squares matrices
   A = [];
   b = [];
   for j = 1:n
      x = p0(1,j);
      y = p0(2,j);
      A = [A; resid_wght(j)*wght(j)*[x 0 y 0 1 0; 0 x 0 y 0 1]];
      b = [b; resid_wght(j)*wght(j)*p1(:,j)];
   end
   
   % Solve for the affine tranform
   aff = pinv(A)*b;
   aff = [aff(1:2) aff(3:4) aff(5:6); 0 0 1];
end

% Determine which constraints were outliers
if nargout >= 2
   pts = aff * [p0; ones(1,n)];
   resid = sqrt(sum((pts(1:2,:) - p1).^2,1))';
   outliers = find(resid > 2.5*sigma);
end

% Determine the mean squared error
if nargout >= 3
   robustError = mean(resid(find(resid <= 2.5*sigma)));
end
