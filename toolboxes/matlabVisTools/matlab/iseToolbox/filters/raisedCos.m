function result = raisedCos(n)
% raisedCos: Returns a raised cosine window of length n
%
result = 0.5 * (1.0 - cos((2*pi/(n-1)) * [0:n-1]'));
