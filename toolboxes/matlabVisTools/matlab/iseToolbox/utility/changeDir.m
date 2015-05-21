function cwd = changeDir(newDir)
%
%AUTHOR: Wandell
%DATE:   12.05.95
%PURPOSE:
%  Duh.
% I got tired of typing this, and the cd command doesn't take a string
%

eval(['cd ',newDir]);

if nargout > 0
 cwd = cd;
end
