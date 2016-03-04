function [ source_eq_image, f_maps_source, f_maps_target,f_maps_source_normalized ] = opp_col_normalization( source_im_name, target_im_name,  rotation_matrix, image_dir, vM_dir, out_dir, opts, varargin)
% color normalization 
% matching the statistics of different channels from source to target
opts_default.matchMethod = 'moments';
opts_default.features.which_features = {'brightness opp','hue opp', 'saturation opp'};
if nargin < 7
    opts = opts_default;
    if nargin < 6
          error(['Need at least 6 inputs: source image name, ', ...
              'target image name, rotation matrix,  img_dir, vonMises dir, and outputdir']);
    end
end 

%matchMethod = optimget(opts,'matchMethod','moments');
matchMethod = opts.matchMethod;
source_im = double(imread(fullfile(image_dir,[source_im_name '.tif'])))./255;
target_im = double(imread(fullfile(image_dir,[target_im_name '.tif'])))./255;
source_stats = load(fullfile(vM_dir,[source_im_name '_stats.mat']));
target_stats = load(fullfile(vM_dir,[target_im_name '_stats.mat']));
opts_matching.source_stats = source_stats.data;
opts_matching.target_stats = target_stats.data;
indx_white_source = logical(source_stats.data.posterior_probs(:,3));
%indx_white_target = target_stats.data.posterior_probs(:,3);

% Color normalization using opponent color space
%% calculate the features in opponent color space
which_features = opts.features.which_features;

% for feature_iter = 1: length(which_features)
%     if (opts.display_progress), fprintf('\nProcessing feature type ''%s'':\n',which_features{feature_iter}); end
%     f_map_source_curr = getFeatures(double(source_im),1,which_features{feature_iter},opts);
%     f_maps_source = cat(3,f_maps_source,f_map_source_curr);
%     f_map_target_curr = getFeatures(double(target_im),1,which_features{feature_iter},opts);
%     f_maps_target = cat(3,f_maps_target,f_map_target_curr);    
% end

f_maps_source = getFeatures(double(source_im),1,which_features,opts);
f_maps_target = getFeatures(double(target_im),1,which_features,opts);
f_maps_source_normalized = cell(1,3);
for feature_iter = 1:length(which_features)
    f_map_source_curr = f_maps_source{feature_iter};
    f_map_target_curr = f_maps_target{feature_iter};
    %% normalization
    if strcmp(matchMethod,'cdf')
        f_map_normalized_curr = matchingCdf(f_map_source_curr, f_map_target_curr, opts_matching);
    elseif strcmp(matchMethod,'moments');
        f_map_normalized_curr = matchingMoments(f_map_source_curr, f_map_target_curr,which_features{feature_iter},opts_matching);
    end 
    f_maps_source_normalized{feature_iter} = f_map_normalized_curr; 
end

%% calculate c2, c3, recover rotated coordinate
% then calculate new rgb
%sat = sqrt(source_rotated(2,:).^2 + source_rotated(3,:).^2); % calculate r
source_rotated_eq = zeros(3,length(f_map_normalized_curr));
source_rotated_eq(1,:) = f_maps_source_normalized{2}; % brightness normalized/equalized
source_rotated_eq(2,:) = f_maps_source_normalized{3}.*cos(f_maps_source_normalized{1}); % c2
source_rotated_eq(3,:) = f_maps_source_normalized{3}.*sin(f_maps_source_normalized{1}); % c3
source_rgb_eq = rotation_matrix.rotation_matrix\source_rotated_eq; % RGB = Rot_matrix \ rotated coordinates
%source_rgb_eq = rotation_matrix_source\source_rotated_eq;
source_rgb_eq_uint8 = uint8(source_rgb_eq*255); 
source_rgb_eq_uint8(:,indx_white_source) = 255;

source_eq_image = reshape(source_rgb_eq_uint8', size(source_im));
% [source_xsize, source_ysize] = size(source_im(:,:,1));
% r = reshape(source_rgb_eq_uint8(1,:),[source_xsize, source_ysize]);
% g = reshape(source_rgb_eq_uint8(2,:),[source_xsize, source_ysize]);
% b = reshape(source_rgb_eq_uint8(3,:),[source_xsize, source_ysize]);
% source_eq_image = cat(3,r,g,b);
%figure; imshow(source_eq_image);
imwrite(source_eq_image,fullfile(out_dir,[source_im_name '.tif']));
%% plot the distribution
% figure;
% subplot(1,3,1);circ_plot(reshape(f_maps_target(:,:,3),[1 numel(r)]),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,3,2);circ_plot(reshape(f_maps_source(:,:,3),[1 numel(r)]),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,3,3);circ_plot(f_maps_target(:,3),'hist',[],40,true,true,'linewidth',2,'color','r');
end

