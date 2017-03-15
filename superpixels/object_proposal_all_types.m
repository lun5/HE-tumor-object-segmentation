function segs = object_proposal_all_types( IMG_DIR, imname, param_string, ...
    max_num_neighbors, num_comps )
% compute the top connected components from the object maps in Burak's
% output. 

% INPUT:
% -IMG_DIR: image and object directory
% -imname: image name
% -param_string: of form _se1_minNuc3_minStr5_minLum9_adjDela, with info
% about minimum nuclear, stroma, and lumen size
% -obj_type: type of object (N-1, S-2, L-3) for finding connected component
% -num_neighbors: number of neighbors to identify the median distances
% between objects (default 15)
% -num_comp: number of top connected component to report (default 10)
% -plot_flag: 1 then plot and save the output (default 0)

% OUTPUT:
% -top_centers: coordinates of centers for top components
% -top_radii: radius of objects for top components
% Luong Nguyen 1/16/2016

if nargin < 5
    num_comps =10;
elseif nargin < 4
    max_num_neighbors = 15;
elseif nargin < 3
    error('Please input: image directory, image name, parameter string');
end
%% LOOP for different thresholds
seg_purple = object_proposal_by_type( IMG_DIR, imname, param_string, ...
    1, max_num_neighbors, num_comps );
seg_white = object_proposal_by_type( IMG_DIR, imname, param_string, ...
    3, max_num_neighbors, num_comps);
% seg_pink = object_proposal_by_type( IMG_DIR, imname, param_string, ...
%     2, max_num_neighbors, num_comps, plot_flag );
% max_purple_cmp = max(seg_purple(:));
% max_pink_cmp = max(seg_pink(:));
% segs = seg_purple + (seg_pink + max_purple_cmp).*(seg_pink > 0) + ...
%     (seg_white + max_purple_cmp + max_pink_cmp + 1).*(seg_white > 0);

% fuse the two masks
max_purple_cmp = max(seg_purple(:));
white_comps = unique(seg_white(:));
segs = seg_purple;
for i = 2:length(white_comps)
    id = white_comps(i);
    mask_white = seg_white == id;
    if sum(seg_purple(mask_white)>0)/sum(mask_white(:)) > 0.5
        mask_diff = ((seg_purple - mask_white) == - 1);
        segs(mask_diff) = max_purple_cmp + id;
    else
        segs(mask_white) = max_purple_cmp + id;
    end
end
% have to do the same for purple
purple_comps = unique(seg_purple(:));
for i = 2:length(purple_comps)
    id = purple_comps(i);
    mask_purple = seg_purple == id;
    if sum(seg_white(mask_purple) >0) /sum(mask_purple(:)) > 0.3
        segs(mask_purple) = id;
    end
end
%% LOOP for different thresholds
% display
%I = imread(fullfile(IMG_DIR,'images',[imname '.tif']));
%bdry = seg2bdry(segs,'imageSize');
% bdry = zeros(size(I,1),size(I,2));
% unique_segs = unique(segs);
% for i = 2:length(unique_segs)
%     id = unique_segs(i);
%     mask = segs == id;
%     bdry = bdry + edge(mask);
% end
% bdry = bdry > 0;
%se = strel('disk',8,4);
%bdry = imdilate(bdry,se);
%bdry_map = double(I)./255.*repmat(~double(bdry),[1 1 3]);
%imwrite(bdry_map,fullfile(outdir,'bdry_im',[imname '.png']));

%segs = segs(1:4:end, 1:4:end); % only with 2kx2k images
segs  = segs + 1;
%parsave(fullfile(outdir,'segmented_images',[imname '.mat']),{segs});

end