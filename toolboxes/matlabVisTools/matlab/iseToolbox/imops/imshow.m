function h=imshow(varargin)
%IMSHOW Display image data.
%   IMSHOW(I) displays the intensity image I.
%
%   IMSHOW(I,[A B]) displays I so that the pixel value A displays
%   as black, the pixel value B displays as white, and pixel
%   values in between display as intermediate shades of gray.
%
%   IMSHOW(RGB) displays the truecolor image RGB.
%
%   IMSHOW(X,MAP) displays the indexed image X,MAP.
%
%   IMSHOW(I,N) displays the intensity image I with a gray scale
%   colormap of length N.  If N is not specified, IMSHOW uses the
%   default value of 256 on a 24-bit display device; otherwise it
%   uses the default value of 64.
%
%   IMSHOW(x,y,A,...) uses the vectors x,y to define the axis 
%   coordinates.
%
%   H = IMSHOW(...) returns a handle to the created image
%   object.
%
%   After displaying the image, if the figure contains only that
%   image and its parent axes, and no other axes or visible
%   uicontrols, then IMSHOW adjusts the axes extent according the
%   current IMBORDER setting and then calls TRUESIZE to adjust
%   the size of the figure and the image. If the current IMBORDER
%   setting is 'tight', then the image will fill the resulting
%   figure. If the current IMBORDER setting is 'loose', then
%   there will be some space between the image axes and the edges
%   of the figure.
%
%   See also IMBORDER, IMAGE, IMAGESC, TRUESIZE, SUBIMAGE, WARP.

%   Clay M. Thompson 5-12-93
%   Revised Steven L. Eddins, April 1996
%   Copyright (c) 1993-1996 by The MathWorks, Inc.
%   $Revision: 5.13 $  $Date: 1996/10/23 19:34:54 $

% 1. Parse input arguments
% 2. Get an axes to plot in.
% 3. Create the image and axes objects and set their display
%    properties.
% 4. If the image is alone in the figure, position the axes
%    according to the current IMBORDER setting and then call
%    TRUESIZE.
%
% Local function:
% ParseInputs
% IsVector
% DoTruesize

[imtype, cdata, cdatamapping, clim, map, xdata, ydata, filename] = ...
    ParseInputs(varargin{:});
imsize = size(cdata);
imsize = imsize(1:2);  % In case ndims(cdata) > 2

axHandle = newplot;
figHandle = get(axHandle, 'Parent');

% Make the image object.
hh = image(xdata, ydata, cdata, ...
        'Parent', axHandle, 'CDataMapping', cdatamapping);

% Set axes and figure properties if necessary to display the 
% image object correctly.
axis image;
set(axHandle, 'XTick', [], 'YTick', [], 'XGrid', 'off', 'YGrid', 'off', ...
        'Visible', 'off');
set(get(axHandle,'Title'),'Visible','on');
set(get(axHandle,'XLabel'),'Visible','on');
set(get(axHandle,'YLabel'),'Visible','on');
if (~isempty(map))
    set(figHandle, 'Colormap', map);
end
if (~isempty(clim))
    set(axHandle, 'CLim', clim);
end

% Do truesize if called for.
if (DoTruesize(figHandle, axHandle))
    if (strcmp(imborder, 'tight'))
        % Have the image fill the figure.
        set(axHandle, 'Units', 'normalized', ...
                'Position', [0 0 1 1]);
    else
        set(axHandle, 'Units', get(figHandle, 'DefaultAxesUnits'), ...
                'Position', get(figHandle, 'DefaultAxesPosition'));
    end
    truesize(figHandle);
end

if (~isempty(filename) & isempty(get(get(axHandle,'Title'),'String')))
    title(filename,'Interpreter','none');
end

if (nargout > 0)
  % Only return handle if caller requested it.
  h = hh;
end

%----------------------------------------------------------------------
% Subfunction ParseInputs
%----------------------------------------------------------------------

function [imtype, cdata, cdatamapping, clim, map, xdata, ...
            ydata, filename] =  ParseInputs(varargin);

filename = '';

if (get(0,'ScreenDepth') > 16)
    defGrayMapLength = 256;
else
    defGrayMapLength = 64;
end

% MAP can be the 2nd or 4th input argument.  If it's empty,
% ignore it.
if ((nargin >= 4) & isempty(varargin{4}))
    warning('Ignoring empty fourth input argument');
    varargin(4) = [];
elseif ((nargin >= 2) & isempty(varargin{2}))
    warning('Ignoring empty second input argument');
    varargin(2) = [];
end

switch length(varargin)
case 0
    error('Not enough input arguments.  See HELP IMSHOW');
    
case 1
    % IMSHOW(I)
    % IMSHOW(RGB)
    % IMSHOW(FILENAME)
    
    if (isstr(varargin{1}))
        % IMSHOW(FILENAME)
        filename = varargin{1};
        [cdata,map] = imread(filename);
        xdata = (1:size(cdata,2));
        ydata = (1:size(cdata,1));
        if (isempty(map))
            if (ndims(cdata) == 3)
                imtype = 'rgb';
            else
                imtype = 'intensity';
            end
            cdatamapping = 'scaled';
            if (isa(cdata, 'double'))
                clim = [0 1];
            elseif (isa(cdata, 'uint8') & ~islogical(cdata))
                clim = [0 255];
            elseif (isa(cdata, 'uint8') & islogical(cdata))
                clim = [0 1];
            end
            map = gray(defGrayMapLength);
            
        else
            imtype = 'indexed';
            cdatamapping = 'direct';
            clim = [];  % irrelevant
        end
    
    elseif (ndims(varargin{1}) == 3)
        % IMSHOW(RGB)
        imtype = 'rgb';
        cdata = varargin{1};
        cdatamapping = 'direct'; % irrelevant for RGB
        clim = [];               % irrelevant for RGB
        map = [];                % irrelevant for RGB
        xdata = [1 size(cdata,2)];
        ydata = [1 size(cdata,1)];
        
    else
        % IMSHOW(I)
        imtype = 'intensity';
        cdata = varargin{1};
        cdatamapping = 'scaled';
        switch class(cdata)
        case 'uint8'
            if (islogical(cdata))
                clim = [0 1];
            else
                clim = [0 255];
            end
        case 'double'
            clim = [0 1];
        otherwise
            error('Unsupported image class');
        end
        map = gray(defGrayMapLength);
        xdata = [1 size(cdata,2)];
        ydata = [1 size(cdata,1)];
        
    end
    
case 2
    % IMSHOW(X,map)
    % IMSHOW(I,N)
    % IMSHOW(I,[a b])
    
    if (prod(size(varargin{2})) == 1)
        % IMSHOW(I,N)
        imtype = 'intensity';
        cdata = varargin{1};
        cdatamapping = 'scaled';
        switch class(cdata)
        case 'uint8'
            if (islogical(cdata))
                clim = [0 1];
            else
                clim = [0 255];
            end
        case 'double'
            clim = [0 1];
        otherwise
            error('Unsupported image class');
        end
        map = gray(varargin{2});
        xdata = [1 size(cdata,2)];
        ydata = [1 size(cdata,1)];
        
    elseif (isequal(size(varargin{2}), [1 2]))
        % IMSHOW(I,[a b])
        imtype = 'intensity';
        cdata = varargin{1};
        cdatamapping = 'scaled';
        clim = varargin{2};
        map = gray(defGrayMapLength);
        xdata = [1 size(cdata,2)];
        ydata = [1 size(cdata,1)];
        
    elseif (size(varargin{2},2) == 3)
        % IMSHOW(X,map)
        imtype = 'indexed';
        cdata = varargin{1};
        cdatamapping = 'direct';
        clim = [];   % irrelevant
        map = varargin{2};
        xdata = [1 size(cdata,2)];
        ydata = [1 size(cdata,1)];
        
    else
        error('Invalid input arguments; see HELP IMSHOW');
        
    end
    
case 3
    % IMSHOW(R,G,B)
    % IMSHOW(x,y,I)
    % IMSHOW(x,y,RGB)
    
    if (ndims(varargin{3}) == 3)
        % IMSHOW(x,y,RGB)
        imtype = 'rgb';
        cdata = varargin{3};
        cdatamapping = 'direct'; % irrelevant
        clim = [];               % irrelevant
        map = [];                % irrelevant
        xdata = varargin{1};
        ydata = varargin{2};
        
    elseif (IsVector(varargin{1}) & IsVector(varargin{2}))
        % IMSHOW(x,y,I)
        imtype = 'intensity';
        cdata = varargin{3};
        cdatamapping = 'scaled';
        switch class(cdata)
        case 'uint8'
            if (islogical(cdata))
                clim = [0 1];
            else
                clim = [0 255];
            end
        case 'double'
            clim = [0 1];
        otherwise
            error('Unsupported image class');
        end
        map = gray(defGrayMapLength);
        xdata = varargin{1};
        ydata = varargin{2};
        
    elseif (isequal(size(varargin{1}), size(varargin{2})) & ...
                isequal(size(varargin{1}), size(varargin{3})))
        % IMSHOW(R,G,B)
        imtype = 'rgb';
        cdata = cat(3,varargin{:});
        cdatamapping = 'direct';        % irrelevant
        clim = [];                      % irrelevant
        map = [];                       % irrelevant
        xdata = [1 size(cdata,2)];
        ydata = [1 size(cdata,1)];
        
    else
        error('Invalid input arguments; see HELP IMSHOW');
        
    end
    
case 4
    % IMSHOW(x,y,X,MAP)
    % IMSHOW(x,y,I,N)
    % IMSHOW(x,y,I,[a b])
    
    if (prod(size(varargin{4})) == 1)
        % IMSHOW(x,y,I,N)
        imtype = 'intensity';
        cdata = varargin{3};
        cdatamapping = 'scaled';
        switch class(cdata)
        case 'uint8'
            if (islogical(cdata))
                clim = [0 1];
            else
                clim = [0 255];
            end
        case 'double'
            clim = [0 1];
        otherwise
            error('Unsupported image class');
        end
        map = gray(varargin{4});
        xdata = varargin{1};
        ydata = varargin{2};
        
    elseif (isequal(size(varargin{4}), [1 2]))
        % IMSHOW(x,y,I,[a b])
        imtype = 'intensity';
        cdata = varargin{3};
        cdatamapping = 'scaled';
        clim = varargin{4};
        map = gray(defGrayMapLength);
        xdata = varargin{1};
        ydata = varargin{2};
        
    elseif (size(varargin{4},2) == 3)
        % IMSHOW(x,y,X,map)
        imtype = 'indexed';
        cdata = varargin{3};
        cdatamapping = 'direct';
        clim = [];                % irrelevant
        map = varargin{4};
        xdata = varargin{1};
        ydata = varargin{2};
        
    else
        error('Invalid input arguments.  See HELP IMSHOW');
        
    end
    
case 5
    % IMSHOW(x,y,R,G,B)
    
    imtype = 'rgb';
    cdata = cat(3,varargin{3:5});
    cdatamapping = 'direct';           % irrelevant
    clim = [];                         % irrelevant
    map = [];                          % irrelevant
    xdata = varargin{1};
    ydata = varargin{2};
    
otherwise
    
    error('Too many input arguments.  See HELP IMSHOW');
    
end

% Catch complex CData case
if (~isreal(cdata))
    warning('Displaying real part of complex input');
    cdata = real(cdata);
end


%%%
%%% Subfunction IsVector
%%%
function tf = IsVector(x)
%ISVECTOR True if x has only one non-singleton dimension.
tf = (sum(size(x)~=1) <= 1);


%%%
%%% Subfunction DoTruesize
%%%
function tf = DoTruesize(figHandle, axHandle)

if (length(findobj(axHandle, 'Type', 'image')) > 1)
    % More than one image in axes
    tf = 0;

else

    figKids = allchild(figHandle);
    kids = [findobj(figKids, 'Type', 'axes') ;
        findobj(figKids, 'Type', 'uicontrol', 'Visible', 'on')];
    if (length(kids) > 1)
        % The axes isn't the only thing around, so don't truesize
        tf = 0;
    else
        % Is axHandle in the default position?
        if (isequal(get(axHandle, 'Position'), ...
                    get(get(axHandle,'Parent'), 'DefaultAxesPosition')))
            % Yes, call truesize
            tf = 1;
            
        else
            % No, don't call truesize
            tf = 0;
        end
    end
end
