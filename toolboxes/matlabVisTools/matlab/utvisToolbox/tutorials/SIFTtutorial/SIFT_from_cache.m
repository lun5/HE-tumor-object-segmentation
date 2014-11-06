function [pos, scale, orient, desc, im, mask] = SIFT_from_cache( im_path, im_name, cache, octaves, intervals )

% [pos, scale, orient, desc, im, mask] = SIFT_from_cache( im_path, im_name, cache, octaves, intervals )
% 
% Apply the SIFT transform to an image specified as a path and and image root
% name, or recover the SIFT transform of the image from cached results, depending
% on the value of the parameter cache.
%
% Input:
% im_path - the path to the image
% im_name - the name of a pgm image, without extension
% cache - if 1 then cached SIFT results are returned if they exist, otherwise
%   they are computed and the cache is recreated.  If 2, then the cache is
%   ignored.  If 3, then the cached versions of Lowe's keypoints are returned.
% octaves - the number of octaves for the SIFT transform (default 4).
% intervals - the number of intervals per ocatve for the SIFT transform (default 2).
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
% im - the image
% mask - the object mask
%
% Thomas F. El-Maraghi
% May 2004

if ~exist('cache')
   cache = 1;
end
if ~exist('octaves')
   octaves = 4;
end
if ~exist('intervals')
   intervals = 2;
end

if exist([im_path,im_name,'.pgm'], 'file')
  im = pgmRead( [im_path,im_name,'.pgm'] )./255;
elseif exist([im_path,im_name,'.ppm'], 'file')
  imRGB = ppmRead( [im_path,im_name,'.ppm'] )./255;
  im = imRGB(:,:,2);
else
  fprintf( 2, ['Cannot find image ' im_name] );
  return;
end

if exist([im_path,im_name,'.mask'])
   mask = pgmRead([im_path,im_name,'.mask'])./255;
else
   mask = ones(size(im));
end

if cache == 1
   if exist([im_path,im_name,'.key.mat'])
      fprintf( 2, 'Reading keypoints for image %s ', im_name );
      [ pos, scale, orient, desc ] = read_keypoint_file( [im_path,im_name] );
      fprintf( 2, '%d keypoints.\n', size(pos,1) );
   else      
      fprintf( 2, 'Computing keypoints for image %s...\n', im_name );
      [ pos, scale, orient, desc ] = SIFT( im, octaves, intervals, mask );   
      write_keypoint_file( [im_path,im_name], pos, scale, orient, desc );
   end
elseif cache == 2
   fprintf( 2, 'Computing keypoints for image %s...\n', im_name );
   [ pos, scale, orient, desc ] = SIFT( im, octaves, intervals, mask );   
   write_keypoint_file( [im_path,im_name], pos, scale, orient, desc );
else      
   fprintf( 2, 'Reading Lowe keypoint file for image %s ', im_name );
   [ pos, scale, orient, desc ] = read_lowe_keypoints( [im_path,im_name], mask );
   fprintf( 2, '%d keypoints.\n', size(pos,1) );
end
