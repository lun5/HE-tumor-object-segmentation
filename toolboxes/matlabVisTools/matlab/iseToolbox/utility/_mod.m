function m = mod(a,n)
%
%  Modulus operator.  
%  Works on matrics, vectors or scalars.

m = a - n .* fix(a./n);
return;

