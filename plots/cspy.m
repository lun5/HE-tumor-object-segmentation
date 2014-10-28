function [] = cspy(varargin)
%CSPY Visualize sparsity pattern.
%   CSPY(S) plots the sparsity pattern of the matrix S with L levels where
%   L by default is number of integers between min matrix value and max
%   matrix value .
%
%   CSPY(S,'Marker', '*') use the this marker to plot the matrix S.
%
%   CSPY(S, 'Marker', {'*', '+'}) use each marker per level in the same
%   order. The size of the marker cell array must be equal to the number of
%   levels, in otherwise is used the first marker for all levels.
%
%   CSPY(S,'MarkerSize', M) use the M markersize to plot the matrix S.
%
%   CSPY(S, 'Levels', L) use L levels to show different colors in the
%   matrix S. If the levels is a vector and has the same size to channels
%   is used one value for channel. Where level hasn't the same size is used
%   the mean value.
%
%   CSPY(S, 'ColorMap', 'hot') use this model color to assign the colors.
%   By default the colormap is 'lines'.
%   
%   CSPY(S, 'XDir', 'Reverse') change the X direction (Reverse or Normal).
%   Default Normal
%
%   CSPY(S, 'YDir', 'Reverse') change the Y direction (Reverse or Normal). 
%   Default Reverse
%   
%   CSPY(S, 'Channels', 'on') use the channel value ('on' or 'off' ) for 
%   plot the images channels in different subplots (rgb for example) or 
%   where the matrix is n x m x c where c are the channels. if the channels
%   value is 'off' and the matrix has c channels with c different to 3(rgb)
%   the channels value is changed to 'on'.
%   
%   Examples:
%   ---------------------------------------------------------
%   cspy(bucky);
%   ---------------------------------------------------------
%   rgb = imread('ngc6543a.jpg');
%   cspy(rgb, 'marker', '.', 'markersize', 5, 'colormap', 'hot'); 
%   ---------------------------------------------------------
%   rgb = imread('ngc6543a.jpg');
%   cspy(rgb, 'marker', '.', 'markersize',5, 'colormap', 'hot', 'levels', 18); 
%   ---------------------------------------------------------
%   rgb = imread('ngc6543a.jpg');
%   cspy(rgb, 'marker', '.', 'markersize',5, 'colormap', 'hot', 'levels', 18, 'channels', 'on'); 
%   ---------------------------------------------------------
%   rgb = imread('ngc6543a.jpg');
%   cspy(rgb, 'marker', '.', 'markersize',5, 'colormap', 'lines', 'levels', [2 15 3], 'channels', 'on');
%   ---------------------------------------------------------
%   rgb = imread('ngc6543a.jpg');
%   cspy(rgb, 'marker', {'.', '+', '*'}, 'markersize',5, 'colormap', 'lines', 'levels', 3, 'channels', 'off');
%
%   Authors:
%   Hugo Gualdron - gualdron@usp.br
%   Jose Rodrigues - junio@icmc.usp.br
%   University of São Paulo

if nargin == 0
    return;
end
islevels = 0;
levels = 0;
MAP = 0;
markerSize = 0;
marker = 0;
ischannels = 0;
xdir = 0;
ydir = 0;

for i=2:2:nargin
    switch(lower(char(varargin(i))))
        case 'marker'
            marker = varargin(i+1);
            if iscell(marker)
                marker = marker{:};
            end
        case 'markersize'
            markerSize = cell2mat(varargin(i+1));
        case 'levels'
            levels = cell2mat(varargin(i+1));
            islevels = 1;
        case 'colormap'
            MAP = char(varargin(i+1));
        case 'channels'
            ischannels = char(varargin(i+1));
        case 'ydir'
            ydir = lower(char(varargin(i+1)));
        case 'xdir'
            xdir = lower(char(varargin(i+1)));
    end
end

matrix = cell2mat(varargin(1));

%% default values
if ~ischannels | ~strcmp('on', ischannels)
    ischannels = 0;
else
    ischannels = 1;
end

channels = size(matrix,3);
if ~ischannels
    %where ischannels is off and the matrix is a rgb image :)
    if channels == 3
        matrix = rgb2gray(matrix);
        channels = 1;
    end
end
plotsdim = ceil(sqrt(channels));
if ~MAP
    MAP = 'lines';
end
if ~markerSize
    markerSize = 10;
end
if isempty(marker) || (~iscell(marker) && ~marker )
    marker = '.';
end

if ~ydir
    ydir = 'reverse';
end

if ~xdir
    xdir = 'normal';
end

if length(levels) > 1 && length(levels) ~= channels
    levels = ceil(mean(levels));
end
%%
for j=1:channels
    if channels > 1
        matrix2d = matrix(:,:,j);
    else
        matrix2d = matrix;
    end
    %matrix2d = flipud(matrix2d);
    %where doesn't elements
    if ~nnz(matrix2d)
        continue;
    end
    if channels > 1
        subplot(plotsdim, plotsdim, j);
    end
    
    [y, x, z] = find(matrix2d);
    minvalue = min(z);
    maxvalue = max(z);
    
    %number of integers between minvalue and maxvalue
    if ~islevels 
        levels = ceil(maxvalue - minvalue +1);
    end
    
    if length(levels) > 1
        clevels = levels(j);
    else
        clevels = levels;
    end
    
    
    step = (maxvalue - minvalue)/clevels;
    colors = eval(strcat(MAP,'(', num2str(clevels), ')'));

    step_init = minvalue;
    step_end = minvalue+step;
    colormap(colors);
    if size(matrix2d,1) == size(matrix2d, 2)
        axis square;
    end
    if clevels > 1 && maxvalue - minvalue >1
        colorbar;
        caxis([minvalue maxvalue]);
    end

    xlim([1, size(matrix2d, 2)]);
    ylim([1, size(matrix2d, 1)]);
    xlabel(['nz = ' num2str(nnz(matrix2d))])
    hold on;
    %t = tic;
    for i=1:clevels
        step_init = minvalue + (i-1)*step;
        step_end = minvalue + i*step;
        if i == clevels
            ids = find(z>=step_init & z<=step_end);
        else
            ids = find(z>=step_init & z<step_end);
        end
        if length(marker) == clevels
            cmarker = char(marker(i));
        else
            cmarker = char(marker(1));
        end
        
        plot(x(ids), y(ids), cmarker, 'MarkerSize', markerSize, 'Color', colors(i,:));
        set(gca,'XDir', xdir);
        set(gca,'YDir', ydir);
        set(gca, 'XAxisLocation', 'top');
        %step_init = step_init + step;
        %step_end = step_end + step;
    end
    %fprintf('cspy time %d seconds\n',toc(t));
end

end
