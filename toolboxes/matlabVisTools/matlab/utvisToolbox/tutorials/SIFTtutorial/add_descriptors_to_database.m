function [database,handle] = add_descriptors_to_database( im, pos, scale, orient, desc, existing_database )

% [database] = add_descriptors_to_database( im, pos, orient, scale, desc, existing_database )
%
% Add an image and its descriptors to a database, or create a new database if
% one is not specified.
% 
% Input:
% im - the image.
% pos, orient, scale, desc - keypoint information from the SIFT function.
% existing_database - database to add the descriptors to.  If omitted, a 
%   new database is created.
%
% Output:
% database - the database structure.
% handle - the integer handle of the descriptor block in the database.
%
% Thomas F. El-Maraghi
% May 2004

if exist('existing_database')
   database = existing_database;
else
   database = empty_descriptor_database;
end

database.num_im = database.num_im + 1;
database.im{database.num_im} = im;
database.index = [database.index; database.num_im*ones(size(pos,1),1)];
if nargout == 2
   handle = database.num_im;
end

database.pos = [database.pos; pos];
database.scale = [database.scale; scale];
database.orient = [database.orient; orient];
database.desc = [database.desc; desc];