function imzoom(varargin)
%IMZOOM Zoom in and out on a 2-D plot.
%   Note: IMZOOM has been grandfathered. IMZOOM calls ZOOM with
%   its input arguments. See HELP ZOOM for usage information.
%
%   See also ZOOM.

%   Clay M. Thompson 1-25-93
%   Revised to call zoom, Steven L. Eddins, September 1996
%   Copyright (c) 1993-1996 by The MathWorks, Inc.
%   $Revision: 5.3 $  $Date: 1996/10/23 19:35:04 $

zoom(varargin{:});
