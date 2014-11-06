function write_keypoint_file( fname, pos, scale, orient, desc )

% write_keypoint_file( fname, pos, orient, scale, desc )
%
% Write a keypoint cache file.
%
% Input:
% fname - the filename of the keypoint file, with path but without the 
%   extension .key.mat
% pos - an Nx2 matrix containing the (x,y) coordinates of the keypoints
%   stored in rows.
% scale - an Nx3 matrix with rows describing the scale of each keypoint 
%   (i.e., first column specifies the octave, second column specifies 
%   the interval, and third column specifies sigma).
% orient - a Nx1 vector containing the orientations of the keypoints 
%   [-pi,pi).
% desc - an Nx128 matrix with rows containing the feature descriptors 
%   corresponding to the keypoints.
%
% Thomas F. El-Maraghi
% May 2004
% Revised June 2005

save([fname,'.key.mat'],'pos','scale','orient','desc');


