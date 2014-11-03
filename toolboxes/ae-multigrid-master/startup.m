% Initialize Matlab environment.
function startup
   % initialize matlab paths
   wd = cd;            % working directory
   addpath(wd);
   % pb
   addpath(fullfile(wd,'grouping','lib'));
   % globalization
   addpath(fullfile(wd,'globalization'));
   % ispc
   addpath(fullfile(wd,'ispc'));
   % grep - added by Luong Nguyen 140929
   addpath('C:\Users\luong_nguyen\Documents\MATLAB\grep04apr06');
end
