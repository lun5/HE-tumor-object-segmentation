%
%	startup file
%

disp([' Running ' ...
      which('startup')]);

clear path;
dirstr='';
% GENERIC STUFF (change these to your directory, unless you are on CS-lab)
global matlabVisRoot;
matlabVisRoot = pwd;%'/h/51/jepson/pub/matlab';

%dirstr=[matlabVisRoot];  % Put the Matlab root on the path
rootstr=[pathsep, matlabVisRoot];  % prefix for subsequent directories

% IMAGES directory (sample images)
%dirstr=[dirstr,rootstr,'/images'];
dirstr = [matlabVisRoot, '/images'];

% ISE TOOLBOX DIRECTORY, TEACHING VERSION
rootstr=[pathsep, matlabVisRoot, '/iseToolbox'];
dirstr=[dirstr,rootstr,'/filters'];
dirstr=[dirstr,rootstr,'/imops'];
dirstr=[dirstr,rootstr,'/matrix'];
dirstr=[dirstr,rootstr,'/pyrTools/MEX'];
dirstr=[dirstr,rootstr,'/pyrTools'];
dirstr=[dirstr,rootstr,'/stats'];
dirstr=[dirstr,rootstr,'/synth'];
dirstr=[dirstr,rootstr,'/tutorials'];
dirstr=[dirstr,rootstr,'/utility'];

% UTVIS TOOLBOX DIRECTORY
rootstr=[pathsep, matlabVisRoot, '/utvisToolbox'];
dirstr=[dirstr,rootstr,'/colour'];
dirstr=[dirstr,rootstr,'/file'];
dirstr=[dirstr,rootstr,'/tutorials'];

path(pathdef);
path([dirstr, pathsep, path]);
clear p dirstr rootstr;
clear startup.m;
