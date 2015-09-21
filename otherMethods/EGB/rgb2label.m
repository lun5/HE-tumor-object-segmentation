function seg = rgb2label(im)
% convert color image back to labels
im_flat = reshape(im,[size(im,1)*size(im,2) 3]);
[~, ~, IC] = unique(im_flat,'rows');
seg = reshape(IC,[size(im,1) size(im,2)]);
end