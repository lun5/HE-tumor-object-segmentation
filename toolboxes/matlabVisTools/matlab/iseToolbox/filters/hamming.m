function result = hamming(n)
% HAMMING: Returns a Hamming window of length n
%
result = 0.54 - 0.46 * cos((2*pi/(n-1)) * [0:n-1]');
