function test3
clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
close all;
format long g;
format compact;
fontSize = 20;

% Read in a standard MATLAB color demo image.
folder = fullfile(matlabroot, '\toolbox\images\imdemos');
baseFileName = 'peppers.png';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
if ~exist(fullFileName, 'file')
	% Didn't find it there.  Check the search path for it.
	fullFileName = baseFileName; % No path this time.
	if ~exist(fullFileName, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end
rgbImage = imread(fullFileName);
% Get the dimensions of the image.  numberOfColorBands should be = 3.
[rows, columns, numberOfColorBands] = size(rgbImage);
% Display the original color image.
subplot(2, 2, 1);
imshow(rgbImage);
title('Original Color Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

% Get hue channel
hsv = rgb2hsv(rgbImage);
h = hsv(:,:,1);
% Display the original color image.
subplot(2, 2, 3);
scaledHue = 255*h;
imshow(scaledHue, []);
title('Hue Image', 'FontSize', fontSize);

% Let's compute and display the histogram.
numberOfBins = 256;
[pixelCount, grayLevels] = hist(h(:), numberOfBins);
subplot(2, 2, 2); 
bar(grayLevels, pixelCount);
grid on;
title('Histogram of Hue Channel', 'FontSize', fontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.

subplot(2, 2, 4); 
cmap = jet(numberOfBins);
r = pixelCount / max(pixelCount(:));
drawWheel(r,cmap)
title('Histogram of Hue Channel', 'FontSize', fontSize);

%=============================================================
function drawWheel(r, cmap)
if (any(r > 1) || any(r < 0))
	error('R must be a vector of values between 0 and 1')
end

if numel(r) ~= size(cmap,1)
	error('Length of r and cmap must be the same')
end

n = numel(r);
innerRadius =  80;
outerRadius = 100;

angles = linspace(0,2*pi,n+1);
newR = innerRadius*(1-r);
% Draw the hue in the annulus.
for k = 1:n
	newR(k);
	%drawSpoke(innerRadius, outerRadius, angles(k), angles(k+1), cmap(k,:));
	drawSpoke(newR(k)    , innerRadius, angles(k), angles(k+1), cmap(k,:));
end

% Draw circle at the center.
line(0,0,'marker','o');
% Draw outer black ring.
%line(cos(angles)*outerRadius, sin(angles)*outerRadius, 'LineWidth', 3, 'Color', 'k');
% Draw inner black ring.
%line(cos(angles)*innerRadius, sin(angles)*innerRadius, 'LineWidth', 3, 'Color', 'k');
axis equal;


%=============================================================
function h = drawSpoke(ri,ro,thetaStart,thetaEnd,c)
xInnerLeft  = cos(thetaStart) * ri;
xInnerRight = cos(thetaEnd)   * ri;
xOuterLeft  = cos(thetaStart) * ro;
xOuterRight = cos(thetaEnd)   * ro;

yInnerLeft  = sin(thetaStart) * ri;
yInnerRight = sin(thetaEnd)   * ri;
yOuterLeft  = sin(thetaStart) * ro;
yOuterRight = sin(thetaEnd)   * ro;

X = [xInnerLeft, xInnerRight, xOuterRight xOuterLeft];
Y = [yInnerLeft, yInnerRight, yOuterRight yOuterLeft];

h = patch(X,Y,c);
set(h,'edgeColor', 'none');

% so that means I can also plot the distributions as lines