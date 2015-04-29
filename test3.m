function test3
clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
close all;
format long g;
format compact;
fontSize = 20;

% % Read in a standard MATLAB color demo image.
% folder = fullfile(matlabroot, '\toolbox\images\imdemos');
% baseFileName = 'peppers.png';
% % Get the full filename, with path prepended.
% fullFileName = fullfile(folder, baseFileName);
% if ~exist(fullFileName, 'file')
% 	% Didn't find it there.  Check the search path for it.
% 	fullFileName = baseFileName; % No path this time.
% 	if ~exist(fullFileName, 'file')
% 		% Still didn't find it.  Alert user.
% 		errorMessage = sprintf('Error: %s does not exist.', fullFileName);
% 		uiwait(warndlg(errorMessage));
% 		return;
% 	end
% end
% rgbImage = imread(fullFileName);
% % Get the dimensions of the image.  numberOfColorBands should be = 3.
% [rows, columns, numberOfColorBands] = size(rgbImage);
% % Display the original color image.
% subplot(2, 2, 1);
% imshow(rgbImage);
% title('Original Color Image', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% 
% % Get hue channel
% hsv = rgb2hsv(rgbImage);
% h = hsv(:,:,1);
% % Display the original color image.
% subplot(2, 2, 3);
% scaledHue = 255*h;
% imshow(scaledHue, []);
% title('Hue Image', 'FontSize', fontSize);

% Let's compute and display the histogram.
mu = 0; kappa = 10;
h = circ_vmrnd(mu,kappa,500);
numberOfBins = 100;%256;
%[pixelCount, grayLevels] = hist(h(:), numberOfBins);
linhist = histogram(h(:), 'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'NumBins',numberOfBins);
pixelCount = linhist.Values;
grayLevels = linhist.BinEdges(2:end) - linhist.BinWidth/2;

angles = grayLevels;
rline = circ_vmpdf(angles,mu,kappa);
hold on;xlim([-pi pi]);
plot(angles,rline,'Color','r','LineStyle','-','LineWidth',3); hold off;
hold off;


%subplot(2, 2, 2); 
% figure;
% bar(grayLevels, pixelCount);
% grid on;
% title('Histogram of Hue Channel', 'FontSize', fontSize);
% %xlim([0 grayLevels(end)]); % Scale x axis manually.
% xlim([-pi pi]);
% hold on;
% % % add the line
% angles = linspace(-pi,pi,numberOfBins);
% rline = circ_vmpdf(angles,mu,kappa);
% plot(angles,rline,'Color','r','LineStyle','-','LineWidth',2); hold off;

%subplot(2, 2, 4); 
h = figure;
%cmap = jet(numberOfBins);
cmap = [0.8 0.8 0.8];
r = pixelCount;% / max(pixelCount(:));
drawWheel(r,linhist.BinEdges,cmap)
hold on
% % add the line
x = cos(angles).*(1 + rline');
y = sin(angles).*(1 + rline');
%plot(x,y,'-k','LineWidth',3); hold off;
PlotAxisAtOrigin(x,y);
title('Histogram of Hue Channel', 'FontSize', fontSize);

%=============================================================
function drawWheel(r, angles, cmap)
% if (any(r > 1) || any(r < 0))
% 	error('R must be a vector of values between 0 and 1')
% end

if size(cmap,1) == 1
    cmap = repmat(cmap,[numel(r),1]);
end

if numel(r) ~= size(cmap,1)
	error('Length of r and cmap must be the same')
end


n = numel(r);
innerRadius =  1;%80;
%outerRadius = 100;

%angles = linspace(-pi,pi,n+1);
newR = innerRadius*(1+r);
% Draw the hue in the annulus.
for k = 1:n
	%newR(k)
	%drawSpoke(innerRadius, outerRadius, angles(k), angles(k+1), cmap(k,:));
	drawSpoke(newR(k)    , innerRadius, angles(k), angles(k+1), cmap(k,:));
end

% Draw circle at the center.
line(0,0,'marker','o');
% Draw outer black ring.
%line(cos(angles)*outerRadius, sin(angles)*outerRadius, 'LineWidth', 3, 'Color', 'k');
% Draw inner black ring.
circle_angles  = linspace(-pi,pi,n+1);
line(cos(circle_angles)*innerRadius, sin(circle_angles)*innerRadius, 'LineWidth', 3, 'Color', 'k');
axis equal;
%axis off;


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