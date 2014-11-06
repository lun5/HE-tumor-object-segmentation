function dist = ...
    mrManifoldDistance(grayM,iSize,numSlices,startPt,[noVal],[radius])
% 
% AUTHOR:  Engel, Wandell
% DATE:    Nov., 1994
% 
% PURPOSE:
%   This code is used to create a mex-file for matlab for the routine
%   mrManifoldDistance().  That routine is used in certain flood
%   fill operations.  It was originally used as part of the
%   unfolding code, but its functionality there has been replaced
%   by mrManDist
%   
%   
%   The input is an array of sample point coordinates that should form
%   a connected manifold in three-space.
%
%   The point of the routine is to compute the distances between a point 
%   in three-space and a set of other points in three-space.  The distance
%   is measured through the connected space of points. 
%
%   DESCRIPTION:
%
%    dist = mrManifoldDistance(grayM,iSize,numSlices,startPt,[noVal],[radius])
%
%   ARGUMENTS:
%    grayM:  A volume of binary data indicating where the gray matter is.
%    iSize:  The size of an image slice through the gray matter
%    numSlices: The number of image slices
%    startPt:   3d coordinates defining where to start the flood fill
%
%   OPTIONAL ARGUMENTS:
%    dimdist:Array of y,x and z separations between points.
%    noVal:  The value returned for unreached locations (default 0)
%    radius: The max distance to flood out
%     (default 0 == flood as far as can)
%
%   RETURNS:
%    dist:  distances to each point from startPt -- same size as grayM
%    nPntsReached:  The number of points reached from the start point.
%   
