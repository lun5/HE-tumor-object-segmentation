function [ source_eq_image, f_maps_source, f_maps_target,f_maps_source_normalized ] = oppColNormalization( source_im, target_im, rotation_matrix, opts, varargin)
opts_default.matchMethod = 'moments';
opts_default.features.which_features = {'hue opp', 'brightness opp','saturation opp'};
if nargin < 4
    opts = opts_default;
    if nargin < 3
          error('Need at least 3 inputs: source image, target image, rotation matrix');
    end
end 

%matchMethod = optimget(opts,'matchMethod','moments');
matchMethod = opts.matchMethod;
% Color normalization using opponent color space
%% calculate the features in opponent color space
which_features = opts.features.which_features;
%rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
%rotation_matrix = rotation_matrix.rotation_matrix;
% for feature_iter = 1: length(which_features)
%     if (opts.display_progress), fprintf('\nProcessing feature type ''%s'':\n',which_features{feature_iter}); end
%     f_map_source_curr = getFeatures(double(source_im),1,which_features{feature_iter},opts);
%     f_maps_source = cat(3,f_maps_source,f_map_source_curr);
%     f_map_target_curr = getFeatures(double(target_im),1,which_features{feature_iter},opts);
%     f_maps_target = cat(3,f_maps_target,f_map_target_curr);    
% end

images = {source_im, target_im};
numClusters = 2; % only purple and pink this time
opts_mixture.noise = 1;
mu_white = 2.24; kappa_white = 30;
f_maps = cell(2,1); save_struct = cell(2,1);
for i = 1:2
   % calculate features
   im_rgb = double(images{i})./255;
   nrows = size(im_rgb,1); ncols = size(im_rgb,2);
   X = reshape(im_rgb,[nrows*ncols,3]);
   rotated_coordinates = rotation_matrix.rotation_matrix*X'; %double([r(:)'; g(:)'; b(:)']);
   mask_white = isolateWhite(im_rgb.*255);
   indx_white = mask_white(:);
   mask_red = isolateRed(im_rgb.*255);
   indx_red = mask_red(:);
   indx_purple_pink = (~indx_white) & (~indx_red);
   
   theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
   sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
   brightness = rotated_coordinates(1,:);
   % von Mises
   X_cart = [cos(theta); sin(theta)]';
   X_cart_pp = X_cart(indx_purple_pink,:); % pp = purple + pink
   %% Call the function
   [ mu_hat_polar_pp,~, kappa_hat_pp,posterior_probs_pp, prior_probs_pp] =...
       moVM_fixWhite(X_cart_pp,numClusters,opts_mixture);
   num_pixels = length(theta);
   posterior_probs = zeros(num_pixels,5);
   posterior_probs(indx_purple_pink,[1:2 5]) = posterior_probs_pp;
   posterior_probs(indx_white,3) = 1;
   posterior_probs(indx_red,4) = 1;
   mu_hat_polar = [mu_hat_polar_pp, mu_white];
   kappa_hat = [kappa_hat_pp, kappa_white];
   prior_probs_pp = prior_probs_pp/num_pixels *sum(indx_purple_pink);
   prior_probs = [prior_probs_pp(1:2), sum(indx_white)/num_pixels,...
       sum(indx_red)/num_pixels, prior_probs_pp(3)];
   save_struct{i} = struct('mu_hat_polar',mu_hat_polar,'kappa_hat',kappa_hat,...
       'posterior_probs',posterior_probs,'prior_probs',prior_probs);
   f_maps{i} = {reshape(theta,[nrows,ncols]), reshape(brightness,[nrows,ncols]),...
       reshape(sat,[nrows,ncols])};
end

opts_matching.source_stats = save_struct{1};
opts_matching.target_stats = save_struct{2};
indx_white_source = logical(opts_matching.source_stats.posterior_probs(:,3));
f_maps_source = f_maps{1}; f_maps_target = f_maps{2};
%f_maps_source = getFeatures(double(source_im),1,which_features,opts);
%f_maps_target = getFeatures(double(target_im),1,which_features,opts);
f_maps_source_normalized = cell(1,3);
for feature_iter = 1:length(which_features)
    f_map_source_curr = f_maps_source{feature_iter};
    f_map_target_curr = f_maps_target{feature_iter};
    %% normalization
    if strcmp(matchMethod,'cdf')
        f_map_normalized_curr = matchingCdf(f_map_source_curr, f_map_target_curr, opts_matching);
    elseif strcmp(matchMethod,'moments');
        f_map_normalized_curr = matchingMoments(f_map_source_curr, f_map_target_curr,which_features{feature_iter}, opts_matching);
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
%source_rgb_eq_uint8(:,indx_white_source) = 255;
source_eq_image = reshape(source_rgb_eq_uint8', size(source_im));

%[source_xsize, source_ysize] = size(source_im(:,:,1));
% r = reshape(source_rgb_eq_uint8(1,:),[source_xsize, source_ysize]);
% g = reshape(source_rgb_eq_uint8(2,:),[source_xsize, source_ysize]);
% b = reshape(source_rgb_eq_uint8(3,:),[source_xsize, source_ysize]);
% source_eq_image = cat(3,r,g,b);
figure; imshow(source_eq_image);

%% plot the distribution
% figure;
% subplot(1,3,1);circ_plot(reshape(f_maps_target(:,:,3),[1 numel(r)]),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,3,2);circ_plot(reshape(f_maps_source(:,:,3),[1 numel(r)]),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,3,3);circ_plot(f_maps_target(:,3),'hist',[],40,true,true,'linewidth',2,'color','r');
end



% [f_maps_source, indx_white_source,indx_red_source] = ...
%     getFeatures(double(source_im)./255,1,which_features,opts);
% [f_maps_target, indx_white_target,indx_red_target] = ...
%     getFeatures(double(target_im)./255,1,which_features,opts);
% indx_purple_pink_source = (~indx_white_source) & (~indx_red_source);
% indx_purple_pink_target = (~indx_white_source) & (~indx_red_source);
% indx_purple_pink = {indx_purple_pink_source,indx_purple_pink_target};
% 
% indx_white = {indx_white_source, indx_white_target};
% indx_red = {indx_red_source, indx_red_target};
% %indx_purple_pink = (~indx_white) && (~indx_red);
% 
% theta = {f_maps_source{1}(:)', f_maps_target{1}(:)'};
% numClusters = 2; % only purple and pink this time
% opts_mixture.noise = 1;
% mu_white = 2.24; kappa_white = 30;
% 
% for i = 1:2 % loop through source and target
%     % 1D von Mises mixture model, freeze, no white version
%     X_cart = [cos(theta{i}); sin(theta{i})]';
%     X_cart_pp = X_cart(indx_purple_pink{i},:); % pp = purple + pink
%     %% Call the function
%     [ mu_hat_polar_pp,~, kappa_hat_pp,posterior_probs_pp, prior_probs_pp] =...
%         moVM_fixWhite(X_cart_pp,numClusters,opts_mixture);
%     num_pixels = length(theta{i});
%     posterior_probs = zeros(num_pixels,5);
%     posterior_probs(indx_purple_pink{i},[1:2 5]) = posterior_probs_pp;
%     posterior_probs(indx_white{i},3) = 1;
%     posterior_probs(indx_red{i},4) = 1;
%     mu_hat_polar = [mu_hat_polar_pp, mu_white];
%     kappa_hat = [kappa_hat_pp, kappa_white];
%     prior_probs_pp = prior_probs_pp/num_pixels *sum(indx_purple_pink);
%     prior_probs = [prior_probs_pp(1:2), sum(indx_white)/num_pixels,...
%         sum(indx_red)/num_pixels, prior_probs_pp(3)];
%     save_struct{i} = struct('mu_hat_polar',mu_hat_polar,'kappa_hat',kappa_hat,...
%        'posterior_probs',posterior_probs,'prior_probs',prior_probs);
% end
