function [ pos, scale, orient, desc ] = read_keypoint_file( fname, mask )

% [ pos, scale, orient, desc ] = read_keypoint_file( fname, mask )
%
% Read a keypoint cache file.
%
% Input:
% fname - the filename of the keypoint file, with path but without the 
%   extension .key.mat
% mask - a mask to filter keypoints before returning them.
%
% Output:
% pos - an Nx2 matrix containing the (x,y) coordinates of the keypoints
%   stored in rows.
% scale - an Nx3 matrix with rows describing the scale of each keypoint (i.e.,
%   first column specifies the octave, second column specifies the interval, and
%   third column specifies sigma).
% orient - a Nx1 vector containing the orientations of the keypoints [-pi,pi).
% desc - an Nx128 matrix with rows containing the feature descriptors 
%   corresponding to the keypoints.
%
% Thomas F. El-Maraghi
% May 2004

load([fname,'.key.mat']);

if exist('mask')
   pts = round(pos);
   keep = zeros(size(pts,1),1);
   for k = 1:size(pts,1)          
      if mask(pts(k,2), pts(k,1)) == 1
         keep(k) = 1;
      end
   end
   k = find(keep == 1);
   pos = pos(k,:);
   scale = scale(k);
   orient = orient(k);
   desc = desc(k,:);
end
