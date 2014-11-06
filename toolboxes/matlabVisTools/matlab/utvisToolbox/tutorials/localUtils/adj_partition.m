% Do least squares fit of filter pair lp, hp
% to constant energy partition (by scaling lp and hp).
function [lp, hp] = adj_partition(lp, hp, N)

 nZero = floor((N+1)/2);
 fS2 = floor( size(lp) ./ 2);
 qw = zeros(N,N);
 qw(nZero+1 - fS2(1): nZero+1+fS2(1), nZero+1 -fS2(2):nZero+1+fS2(2)) = lp;
 alp = abs(fftshift(fft2(fftshift(qw)))).^2;

 qw = zeros(N,N);
 qw(nZero+1 - fS2(1): nZero+1+fS2(1), nZero+1 -fS2(2):nZero+1+fS2(2)) = hp;
 ahp = abs(fftshift(fft2(fftshift(qw)))).^2;

 A = zeros(2, 2);
 b = zeros( 2, 1);
 A(1,1) = sum(sum(alp.^2));
 A(2,1) = sum(sum(alp .* ahp));
 A(1,2) = A(2,1);
 A(2,2) = sum(sum(ahp.^2));
 b( 1,1) = sum(sum(alp));
 b( 2,1) = sum(sum(ahp));
 x = A\b;
 x = x .^ 0.5;
 fprintf(1, 'Scaled lp and hp by %f %f\n', x(1), x(2));
 lp = x(1) * lp;
 hp = x(2) * hp;
