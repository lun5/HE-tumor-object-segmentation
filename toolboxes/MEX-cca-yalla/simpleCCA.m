A = rand(10,10);
% make the matrix symmetric
A = A*A';
% make the values go between 0 and 1
A = A/max(A(:));
%A = A + spdiags(ones(size(A,1)),0,size(A,1),size(A,2));
%A = A.*A>0.5;

% A = should be a sparse affinity matrix, otherwise force it to be sparse
% as done below. tau = threshold
tau = 0.7;
[lblA zA] = cca(sparse(A),tau,1);
lbl = lblA;
z = zA;

% another way to reproduce the same result
B = A > tau;
[lblB zB] = cca(sparse(B),0,1);

% confirm the results are the same
% the index of numbers tells us the cluster # for each node
[lblA lblB]


% some more play with the variables
% each pixel is given a label stored in: lbl
ulbl = unique(lbl);
% total number of connected components
% the output variable z carries the same information
length(ulbl)

% sort the label information
ulblLen = [];
for t = 1:length(ulbl)
  ulblLen(t) = length(find(lbl == ulbl(t)));
end
[sulblLen,iulblLen] = sort(-ulblLen);
sulblLen = -sulblLen;
