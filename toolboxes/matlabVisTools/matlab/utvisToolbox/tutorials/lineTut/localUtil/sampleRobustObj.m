function [E, tBins, rBins] = sampleRobustObj(pts, x0, sigma, ntBin, nrBin, ...
                                             rRange, nonMaxRadius)
% sampleRobustObj: 
%
% [E tBins rBins] = sampleRobustObj(pts, x0, sigma, ntBin, nrBin, rRange,...
%                                   nonMaxSupp)

x0 = x0(:);  % ensure it is a column vector

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

% Compute bins in radius
rBins = (0:(nrBin-1))'*(rRange(2)-rRange(1))/nrBin + rRange(1);
dr = rBins(2)-rBins(1);

% Compute normals given theta
nt = [cos(tBins) sin(tBins)];

% Allocate objective function array.
E = zeros(length(rBins), length(tBins));
rE = repmat(rBins, 1, length(tBins));

% Add contribution to E for each edgel
for k = 1:size(pts,1)
  
  % Solve n(theta) dot pt + r = 0 for r:
  r = - nt * (pts(k,:)' - x0);
  
  err = rE - repmat(r', size(rE,1), 1);
  err = err.^2;
  E = E+ err ./(sigma^2 + err);
end
  
if nonMaxRadius>0
  E = -localMax(-E, nonMaxRadius, -size(pts,1));
end

%% Crop out extra columns from E used for non-max suppression
if rad > 0
  E = E(:, (rad+1):(end-rad));
  tBins = tBins((rad+1):(end-rad));
end

return;
