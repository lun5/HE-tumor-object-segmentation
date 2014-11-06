function lplot(x)
% lplot(x): displays a vector as a discrete-time signal

index = 0:length(x)-1;

msize = size(x);
if( msize(1) > msize(2))
  x = x';
end

N = length(x);
xx = [zeros(1,N);x;zeros(1,N)];
indexis = [index;index;index];
xdiscrete = [0 xx(:)' 0];
idiscrete = [-1 indexis(:)' N];

plot(idiscrete, xdiscrete, index, x, 'o');
return
