function [IM_IDX, TRANS, THETA, RHO, DESC_IDX, NN_IDX, WGHT] = hough_similarity( database, pos, scale, orient, desc, threshold, do_nonmax )

% [im_idx trans rot rho idesc inn wght] = hough_similarity( database, pos, orient, scale, desc, threshold, do_nonmax )
%
% Compute a similarity transform between a test set of descriptors and an
% image in the database using the hough transform.  
%
% This function is based on portions of:
% [1] David G. Lowe, "Local feature view clustering for 3D object recognition", 
%     IEEE Conference on Computer Vision and Pattern Recognition, Kauai, Hawaii 
%     (December 2001), pp. 682-688. 
%
% Input:
% database - the descriptor database.
% pos, orient, scale, desc - output of the SIFT function for the image
%   to match to the database.
% threshold - minimum number of matches before a model is selected.
% do_nonmax - true if nonmax suppression should be performed on the
%   histogram before returning candidate transforms.
%
% Output:
% im_idx - index of the image in the database that matches the test image.
% trans - tranlation from the test image to the matching database image.
% rot - rotation from the test image to the matching database image.
% rho - scaling from the test image to the matching database image.
% idesc - indices of the descriptors in the test image that agree
%   with the similarity transform.
% inn - indices of the nearest neighbours corresponding to the keypoints
%   specified by idesc.
% wght - the weight of the hough histogram bin for the transform.
%
% Thomas F. El-Maraghi
% May 2004

if ~exist('threshold')
   threshold = 1.5;
end
if ~exist('do_nonmax')
   do_nonmax = 0;
end

% Get the nearest neigbours from the database.
nn = find_nearest_neighbours( database, desc );
desc_idx = find(nn>0);
pos = pos(desc_idx,:);
scale = scale(desc_idx);
orient = orient(desc_idx);
n = nn(desc_idx);
nn_index = database.index(n);
nn_pos = database.pos(n,:);
nn_scale = database.scale(n);
nn_orient = database.orient(n);

% Set up the orientation bins
rot_bin_spacing = 2*pi/12;
rot_bins = -pi:rot_bin_spacing:(pi - rot_bin_spacing);
num_rot_bins = length(rot_bins);

% Set up the scale bins
scale_bins = 2.^[-8:1:8];
num_scale_bins = length(scale_bins);

IM_IDX = [];
TRANS = [];
THETA = [];
RHO = [];
DESC_IDX = cell(1);
NN_IDX = cell(1);
WGHT = [];
n_transforms = 0;

% Loop over all of the models in the database.
for i_image = 1:database.num_im
   model_matches = find(nn_index == i_image)';
   if length(model_matches) > threshold
      
      % Set up the position bins
      trans_bin_spacing = 0.25*max(size(database.im{i_image}));
      X = [-10*trans_bin_spacing:trans_bin_spacing:10*trans_bin_spacing];
      [x_coords y_coords] = meshgrid( X );
      trans_bins = [x_coords(:) y_coords(:)]';
      num_trans_bins = length(X);      
      
      % Initialize the histogram to zero
		hough_hist = zeros(num_trans_bins,num_trans_bins,num_rot_bins,num_scale_bins);
      keypoints = cell(size(hough_hist));
      nn_idx = cell(size(hough_hist));
      
      for k = model_matches
         x1 = pos(k,:);
         x2 = x1 + 10*scale(k)*[cos(orient(k)) sin(orient(k))];
         
         n_x1 = nn_pos(k,:);
         n_x2 = n_x1 + 10*nn_scale(k)*[cos(nn_orient(k)) sin(nn_orient(k))];
         
         A = [ x1(1) -x1(2) 1 0; ...
               x1(2)  x1(1) 0 1; ...
               x2(1) -x2(2) 1 0; ...
               x2(2)  x2(1) 0 1 ];
         b = [ n_x1'; n_x2' ];
         
         c = pinv(A)*b;
         
         theta = atan2( c(2), c(1) );
         rho = c(1) / cos(theta);
         trans = c(3:4);
         
         % Histogram the translation into the four nearest bins
         trans_diff = sqrt(sum((trans_bins - repmat(trans,1,num_trans_bins*num_trans_bins)).^2,1));   
         for k_trans = 1:4
            [trans_val i_trans_bin] = min(trans_diff);
            trans_diff(i_trans_bin) = max(trans_diff);
            trans_wght = prod( max(1 - abs(trans-trans_bins(:,i_trans_bin))/trans_bin_spacing,0) );
            
            % Histogram the rotation into the 2 nearest bins
            rot_diff = abs(mod( theta - rot_bins + pi, 2*pi ) - pi);
            for k_rot = 1:2
               [rot_val i_rot_bin] = min(rot_diff);
               rot_diff(i_rot_bin) = max(rot_diff);
               rot_wght = max(1 - rot_val/(1.5*rot_bin_spacing),0);
               
               vote_wght = trans_wght * rot_wght;
               if vote_wght > 0.0 
                  % Histogram the scaling into the 2 nearest bins
                  rho_diff = abs(rho - scale_bins);
                  for k_scale = 1:2
                     [scale_val i_scale_bin] = min(rho_diff);
                     rho_diff(i_scale_bin) = max(rho_diff);
                     
                     % Split the x and y indices of the tranlation bin
                     i_trans = get_multiple_indices( i_trans_bin, [num_trans_bins num_trans_bins] );
                     
                     % Accumulate the histogram and record which points were placed in which bin
                     hough_hist(i_trans(1),i_trans(2),i_rot_bin,i_scale_bin) = hough_hist(i_trans(1),i_trans(2),i_rot_bin,i_scale_bin) + vote_wght;
                     keypoints{i_trans(1),i_trans(2),i_rot_bin,i_scale_bin} = [keypoints{i_trans(1),i_trans(2),i_rot_bin,i_scale_bin}; desc_idx(k)];
                     nn_idx{i_trans(1),i_trans(2),i_rot_bin,i_scale_bin} = [nn_idx{i_trans(1),i_trans(2),i_rot_bin,i_scale_bin}; n(k)];
                  end
               end
            end            
         end         
      end
      
      % Perform nonmax suppression on the histogram
      if do_nonmax
         h = hough_hist;
         for ix = 1:num_trans_bins
            for iy = 1:num_trans_bins
               for r = 1:num_rot_bins
                  for s = 1:num_scale_bins
                     val = h(ix,iy,r,s);
                     if val > threshold
                        is_max = 1;
                        for ixs = max(min(ix+[-1 0 1],num_trans_bins),1)
                           for iys = max(min(iy+[-1 0 1],num_trans_bins),1)
                              for rs = mod(r+[-1 0 1]+num_rot_bins-1,num_rot_bins)+1
                                 for ss = max(min(s+[-1 0 1],num_scale_bins),1)
                                    if val < h(ixs,iys,rs,ss)
                                       is_max = 0;
                                    end
                                 end
                              end
                           end
                        end
                        if is_max == 0
                           hough_hist(ix,iy,r,s) = 0;
                        end
                     else
                        hough_hist(ix,iy,r,s) = 0;
                     end
                  end
               end
            end
         end                         
      end
      
      % Return any candidate transformations
      candidate_transforms = find(hough_hist(:) > threshold)';
      for candidate = candidate_transforms
         n_transforms = n_transforms + 1;
         
         idx = get_multiple_indices( candidate, [num_trans_bins num_trans_bins num_rot_bins num_scale_bins] );
         
         IM_IDX = [IM_IDX; i_image];
         TRANS = [TRANS; trans_bins(:,get_linear_index(idx(1:2),[num_trans_bins num_trans_bins]))'];
         THETA = [THETA; rot_bins(idx(3))];
         RHO = [RHO; scale_bins(idx(4))];
         DESC_IDX{n_transforms} = keypoints{idx(1),idx(2),idx(3),idx(4)};
         NN_IDX{n_transforms} = nn_idx{idx(1),idx(2),idx(3),idx(4)};
         WGHT = [WGHT; hough_hist(idx(1),idx(2),idx(3),idx(4))];
      end                     
   end
end
