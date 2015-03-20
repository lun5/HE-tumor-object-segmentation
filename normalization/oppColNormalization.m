% current trouble: I have not figured out how to match circular variables
function [ source_eq_image ] = oppColNormalization( source_im, target_im, rotation_matrix, opts, varargin)

if nargin < 4
    opts = [];
    if nargin < 3
          error('Need at least 3 inputs: source image, target image, rotation matrix');
    end
end 

%matchMethod = optimget(opts,'matchMethod','moments');
matchMethod = opts.matchMethod;
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
f_maps_source_normalized = zeros(length(which_features),numel(source_im(:,:,1)));
for feature_iter = 1:length(which_features)
    f_map_source_curr = f_maps_source{feature_iter};
    f_map_target_curr = f_maps_target{feature_iter};
    %% normalization
    if strcmp(matchMethod,'cdf')
        f_map_normalized_curr = matchingCdf(f_map_source_curr, f_map_target_curr);
    elseif strcmp(matchMethod,'moments');
        f_map_normalized_curr = matchingMoments(f_map_source_curr, f_map_target_curr,which_features{feature_iter});
    end 
    f_maps_source_normalized(feature_iter,:) = f_map_normalized_curr; 
end

%% calculate c2, c3, recover rotated coordinate
% then calculate new rgb
%sat = sqrt(source_rotated(2,:).^2 + source_rotated(3,:).^2); % calculate r
source_rotated_eq = zeros(3,length(f_map_normalized_curr));
source_rotated_eq(1,:) = f_maps_source_normalized(2,:); % brightness normalized/equalized
source_rotated_eq(2,:) = f_maps_source_normalized(3,:).*cos(f_maps_source_normalized(1,:)); % c2
source_rotated_eq(3,:) = f_maps_source_normalized(3,:).*sin(f_maps_source_normalized(1,:)); % c3
source_rgb_eq = rotation_matrix.rotation_matrix\source_rotated_eq; % RGB = Rot_matrix \ rotated coordinates
%source_rgb_eq = rotation_matrix_source\source_rotated_eq;
source_rgb_eq_uint8 = uint8(source_rgb_eq*255); 
[source_xsize, source_ysize] = size(source_im(:,:,1));
r = reshape(source_rgb_eq_uint8(1,:),[source_xsize, source_ysize]);
g = reshape(source_rgb_eq_uint8(2,:),[source_xsize, source_ysize]);
b = reshape(source_rgb_eq_uint8(3,:),[source_xsize, source_ysize]);

source_eq_image = cat(3,r,g,b);
figure; imshow(source_eq_image);

%% plot the distribution
% figure;
% subplot(1,3,1);circ_plot(reshape(f_maps_target(:,:,3),[1 numel(r)]),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,3,2);circ_plot(reshape(f_maps_source(:,:,3),[1 numel(r)]),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,3,3);circ_plot(f_maps_target(:,3),'hist',[],40,true,true,'linewidth',2,'color','r');
end

