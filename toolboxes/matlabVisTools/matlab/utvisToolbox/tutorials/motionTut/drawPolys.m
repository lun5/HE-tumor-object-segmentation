function [] = drawPolys(fig, points, penColour, im, sclPix)
%% [] = drawPolys(fig, points, penColour, im, sclPix)
%% Draw polygonal path in figure with handle fig, through vertices:
%%   points = [x(1) x(2) ... x(n); y(1) y(2) ... y(n)];
%% Use x(k) = nan to indicate a break between successive polygonal
%% paths.
%% OPTIONAL: 
%%   penColour = [r g b] (default [1 0 0] (red))
%%   im   Clear figure, draw image in figure

   if nargin < 5
     sclPix = 1; 
   end
   if nargin >= 4
     figure(fig); close; figure(fig);
     %% Display image
     sz = size(im);
     image(im);
     if (size(sz,2) == 2)
       colormap(gray(256)); 
     end
     resizeImageFig(fig, sz, sclPix);
   else
     figure(fig);
   end
   if nargin < 3
     penColour = [1 0 0];
   end

   hold on;
   prev = [];
   for k=size(points,2):-1:1
     current = points(:,k);
     if isnan(current(1))
       current = [];
     end

     if size(prev,1) > 0 & size(current,1) > 0 
       % Plot k-th line from x0,y0 to x1,y1 
       plot([prev(1) current(1)], [prev(2) current(2)], 'LineWidth', 1, ...
            'Color', penColour);
     end
     prev = current;

   end

   return;