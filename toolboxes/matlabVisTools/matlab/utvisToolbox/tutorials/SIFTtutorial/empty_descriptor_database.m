function [database] = create_empty_descriptor_database()

% [database] = create_empty_descriptor_database()
%
% Create an empty descriptor database.
%
% Thomas F. El-Maraghi
% May 2004

database.pos = [];
database.scale = [];
database.orient = [];
database.desc = [];
database.index = [];
database.im = {};
database.num_im = 0;
