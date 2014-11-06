function ffilt = fcos_wave(N)
% Filters L(k), H(k), L0(k), H0(k) for self-inverting
% multiscale pyramid. (CSC2503 Lecture 4.98)
%
% The frequency 0 is returned in pixel position nZero = floor(N/2)+1;
%
% ffilt(1,:)= flp = FT of low-pass filter = L(k),
% ffilt(2,:)= fhp = FT of high-pass filter = H(k),
% ffilt(3,:)= flp0 = FT of initial low-pass filter = L0(k),
% ffilt(4,:)= fhp0 = FT of initial high-pass filter = H0(k).

 ffilt = zeros(4,N);
 
 nZero = floor(N/2)+1;
 nMax = N - nZero;
 nMin = -nZero + 1;
 nInt = floor(N/4);
 nLow = floor(N/8);
 
 x = [nMin:1:nMax];
 x(nZero:N) = min(max(x(nZero:N) - nInt, 0), nInt);
 x(1:nZero) = max(min(x(1:nZero) + nInt, 0), -nInt);
 x= (pi/(2*nInt)) * x;

 ffilt(3,:) = cos(x);
 ffilt(4,:) = sin(abs(x));
 
 x = [nMin:1:nMax];
 x(nZero:N) = min(max(x(nZero:N) - nLow, 0), nLow);
 x(1:nZero) = max(min(x(1:nZero) + nLow, 0), -nLow);
 ffilt(1,:) = cos((pi/(2*nLow)) * x);
 ffilt(2,:) = sin((pi/(2*nLow)) * abs(x));

 
 
