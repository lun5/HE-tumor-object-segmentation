% testing the speed of affinity calculation
% feature type = luminance; affinity type = difference
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
opts_input = setEnvironment_inputs;

I = imread(fullfile(tiles_dir, '036YaafqbAvH5su.tif'));
imshow(I); 

smallestSize = 128;
rows = size(I,1); cols = size(I,2);
% take blocks that are 2, 4,8, 16 times of this block.
opts_affinity = setEnvironment_affinity;
for mult = [1,2,4,6,8,10,12,14,16]
    blockSizeR = smallestSize*mult;
    blockSizeC = smallestSize*mult;
    wholeBlockR = floor(rows/ blockSizeR);
    wholeBlockC = floor(cols/ blockSizeR);
    row1 = 1;
    row2 = row1 + blockSizeR - 1;
    col1 = 1;
    col2 = col1 + blockSizeC - 1;
    cropped_im = I(row1:row2, col1:col2,:);
    fprintf('Size of the block is %d\n',blockSizeR);
    [affinity_matrix, im_sizes] = calculateAffinity(cropped_im, opts_affinity);
end
