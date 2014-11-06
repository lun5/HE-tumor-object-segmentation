function [index] = get_linear_index( indices, dims )

% [index] = get_linear_index( indices, dims )
%
% Compute the linear index given the multidimensional index.
%
% Input:
% indices - the multidimensional indices.
% dim - vector of dimensions of the matrix.
%
% Output:
% index - the linear index.
%
% Thomas F. El-Maraghi
% May 2004

index = indices(1);
for k = 2:length(indices)
   index = index + prod(dims(1:(k-1)))*(indices(k)-1);
end
