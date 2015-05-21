function [indices] = get_multiple_indices( index, dims )

% [indices] = get_multiple_indices( index, dims )
%
% Compute the multidimensional indices from a linear index.
%
% Input:
% index - the linear index.
% dim - vector of dimensions of the matrix.
%
% Output:
% indices - the multidimensional indices.
%
% Thomas F. El-Maraghi
% May 2004

index = index - 1;
for k = length(dims):-1:1
   if k > 1
      divisor = prod(dims(1:(k-1)));
      tmp = floor( index / divisor );
      index = index - divisor*tmp;
      indices(k) = tmp + 1;
   else
      indices(k) = index + 1;
   end
end
