function [ h ] = place_image( image, location, param )
%places an image centered at the middle most pixel at the coordinate
%specified by the 1d vector location
%
%by default, all zero pixels are set to tranparent
%
%returns the handle for the placed image

%grj 10/22/14

if ~exist('param', 'var')
    param = [];
end

param = ml_initparam(param, struct( ...
        'subsize', 400, ...
        'alpha', [0,0,0] ...
        ));

%scale such maximum is 1
imsize = size(image);
% image = image./max(image(:));

image(image > 1) = 1;
image(image < 0) = 0;

h = imagesc(((1:1:imsize(1)) - imsize(1)/2)/param.subsize + location(1),((1:1:imsize(2))-imsize(2)/2)/param.subsize + location(2), image);

alphamap = zeros(size(image));
for i = 1:length(param.alpha)
    alphamap(:,:,i) = param.alpha(i);
end

set(h, 'AlphaDataMapping', 'none')
set(h, 'AlphaData', ~all(image == alphamap,3))


end