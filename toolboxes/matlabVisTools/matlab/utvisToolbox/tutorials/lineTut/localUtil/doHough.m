function [H, tBins, rBins] = doHough(pts, x0, ntBin, nrBin, rRange, nonMaxRadius)
% [H tBins rBins] = doHough(p, x0, ntBin, nrBin, rRange, nonMaxRadius)

% ensure it is a column vector
x0 = x0(:);  

% Compute bins in theta, in [0,pi). 
rad = nonMaxRadius;
% If rad = nonMaxRadius > 0, add rad extra bins on each
% end of the theta range for use during non-max suppression.
tBins = (0:(ntBin-1))'*pi/ntBin;  % tBins sets the lower side of bin.
dt = tBins(2)-tBins(1);
if rad > 0
  tBins = [tBins(1)+dt*(-rad:-1)'; tBins]; % Add extra bins below theta=0
  tBins = [tBins; tBins(end)+dt*(1:rad)']; % Add extra bins beyond theta=pi
end

% Compute bins for r in [rRange(1) rRange(2)).
% Add one extra bin at end of rBins in r for out-of-range.
rBins = (0:nrBin)'*(rRange(2)-rRange(1))/nrBin + rRange(1);
dr = rBins(2)-rBins(1);

% Compute normals given theta
nt = [cos(tBins) sin(tBins)];

% Allocate Hough array.
H = zeros(length(rBins), length(tBins));
for k = 1:size(pts,1)
  
  % Solve n(theta) dot pt + r = 0 for r:
  r = - nt * (pts(k,:)' - x0);
  
  % Quantize r
  qr = floor((r - rBins(1))/dr) + 1;
  % Crop out-of-range values of r, give them the index = length(rBins)
  qr(qr>length(rBins)) = length(rBins);
  qr(qr<=0) = length(rBins);
  
  % Increment Hough transform
  idx = (0:(length(tBins)-1))'* length(rBins) + qr;
  H(idx) = H(idx)+1;
  
end

%% Crop out extra row used for indicating r is out-of-range
H = H(1:end-1, :);
rBins = rBins(1:end-1);

%% Do non-max suppression, if desired
if nonMaxRadius > 0
  H = localMax(H, nonMaxRadius);
end

%% Crop out extra columns from H used for non-max suppression
if rad > 0
  H = H(:, (rad+1):(end-rad));
  tBins = tBins((rad+1):(end-rad));
end

return;

