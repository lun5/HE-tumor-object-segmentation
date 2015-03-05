% imStats(IM1,IM2)
%
% Report image (matrix) statistics.
% When called on a single image IM1, report min, max, mean, stdev, 
% and kurtosis.
% When called on two images (IM1 and IM2), report min, max, mean, 
% stdev of the difference, and also SNR (relative to IM1).

% Eero Simoncelli, 6/96.

function [] = imStats(im1,im2)

if (exist('im2') == 1)
  difference = im1 - im2;
  if isreal(difference)
    [mn,mx] = range2(difference);
  else
    [rmn,rmx] = range2(real(difference));
    [imn,imx] = range2(imag(difference));
    mn = rmn+i*imn; mx = rmx+i*imx;
  end
  mean = mean2(difference);
  v = var2(difference,mean);
  if (v < realmin) 
    snr = Inf;
  else
    snr = 10 * log10(var2(im1)/v);
  end
  fprintf(1, 'Difference statistics:\n');
  fprintf(1, '  Range: [%c, %c]\n',mn,mx);
  fprintf(1, '  Mean: %f,  Stdev (rmse): %f,  SNR (dB): %f\n',...
      mean,sqrt(v),snr);
else
  if isreal(im1)
    [mn,mx] = range2(im1);
  else
    [rmn,rmx] = range2(real(im1));
    [imn,imx] = range2(imag(im1));
    mn = rmn+i*imn; mx = rmx+i*imx;
  end
  mean = mean2(im1);
  stdev = sqrt(var2(im1,mean));
  kurt = kurt2(im1, mean, stdev^2);
  fprintf(1, 'Image statistics:\n');
  fprintf(1, '  Range: [%f, %f]\n',mn,mx);
  fprintf(1, '  Mean: %f,  Stdev: %f,  Kurtosis: %f\n',mean,stdev,kurt);
end
  
