%% proxy for contours2ucm since I do not want output display
% Luong Nguyen
function [ucm] = proxy_contours2ucm(gPb_orient,fmt)
if nargin<2, fmt = 'doubleSize'; end
[~, ucm] = evalc(['contours2ucm(gPb_orient,''' fmt ''');']);
end