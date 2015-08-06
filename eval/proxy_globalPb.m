%% proxy for globalPb since I do not want output display
% Luong Nguyen
function [gPb_orient] = proxy_globalPb(imgFile, outFile, rsz)
if nargin<3, rsz = 1.0; end
[~, gPb_orient] = evalc(['globalPb(''' imgFile ''',''' outFile ''',' num2str(rsz) ');' ]);
end
