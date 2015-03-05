% V = VAR2(MTX,MEAN)
%
% Sample variance of a matrix.
%  MEAN (optional) makes the calculation faster.

function res = var2(mtx, mn)

if (exist('mn') ~= 1)
	mn =  mean(mean(mtx));
end
	
res = sum(sum(abs(mtx-mn).^2)) / (prod(size(mtx)) - 1);
