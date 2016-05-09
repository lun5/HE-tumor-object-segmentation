function [top_centers, top_radii] = top_connected_comp( IMG_DIR, outdir, imname, param_string, ...
    obj_type, max_num_neighbors, num_comps, plot_flag )
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

if nargin < 7
    plot_flag = 0;
elseif nargin < 6
    num_comps =10;
elseif nargin < 5
    max_num_neighbors = 15;
elseif nargin < 4
    obj_type = 1;
elseif nargin < 3
    error('Please input: image directory, image name, parameter string');
end

% read in image, mask, adjacency map
%I = double(imread(fullfile(IMG_DIR,'images', [imname '.jpg'])));
%nrow = size(I,1); ncol = size(I,2);
info = imfinfo(fullfile(IMG_DIR,'images',[imname '.tif']));
nrow = info.Height; ncol = info.Width; 
adj_map = dlmread(fullfile(IMG_DIR,'adjDela_new',[imname param_string '_adjDela']),',',0,0);
%circle_map = dlmread(fullfile(IMG_DIR,'circle_map',[imname param_string '_circle_map']),',',1,0);
num_total_objects = adj_map(1,1);
adj_map = adj_map(2:end,:);
obj_types = adj_map(:,3);
obj_coords = adj_map(:,[5 4]);
obj_radii = sqrt(adj_map(:,2)./pi);

% find the centers and radii of the object of certain type
indx_type = obj_types == obj_type;
% centers = obj_coords(indx_type,:);
% radii = obj_radii(indx_type);
num_superpx = sum(indx_type);

adj_obj = adj_map(indx_type,:); 
adj_obj_coords = adj_obj(:,[5 4]);
adj_obj_radii = sqrt(adj_obj(:,2)./pi);
%adj_matrix = sparse(num_superpx, num_superpx);
%adj_dist = sparse(num_superpx, num_superpx);
% figure; imshow(I); hold on;
% for i = 1:size(adj_obj)
%     viscircles(adj_obj_coords(i,:),adj_obj_radii(i));
% end
% hold off;
% slow way can be fast if we calculate the index first
indx_1 = cell(num_superpx,1);
indx_2 = cell(num_superpx,1);
dist_values = cell(num_superpx,1);
%proportion_features = cell(num_superpx,1); 
%tic;
first_indx = find(obj_types == obj_type,1,'first');
for i = 1:num_superpx
   id1 = adj_obj(i,1);
   num_neighbors_input = adj_obj(i,6);% number of neighbors
   %proportion_features{i} = [0 0];
   if num_neighbors_input > 0
       % neighbors of neighbors, neighbor-ception
       n_n = adj_map(adj_obj(i,7:(6+num_neighbors_input)),[1,7:end]);
       n_n = n_n(:)'; n_n(n_n == 0) = []; 
       neighbor_indx = [adj_obj(i,7:(6+num_neighbors_input)), n_n];
       neighbor_indx(neighbor_indx == id1) = [];
       %neighbor_indx = adj_obj(i,7:(6+num_neighbors_input));
       
       % check if the neighbor is actually nuclei/whatever
       neighbor_types = obj_types(neighbor_indx);
       % calculate feature vector: proportion of purple + pink
       %proportion_features{i} = [sum(neighbor_types == 1) sum(neighbor_types == 2)]./length(neighbor_types);
       neighbors_of_same_type = neighbor_indx(neighbor_types == obj_type);
       num_neighbors_of_same_type = length(neighbors_of_same_type);
       % limit the number of neighbor to be fewer than max_num_neighbors
       %num_neighbors_of_same_type = min(max_num_neighbors, num_neighbors_of_same_type);
       if num_neighbors_of_same_type > 0
           %indx_1{i} = repmat(id1,[1, num_neighbors_of_same_type]);
           %indx_2{i} = neighbors_of_same_type(1:num_neighbors_of_same_type);
           distances = zeros(1, num_neighbors_of_same_type);
           for j = 1:num_neighbors_of_same_type
               distances(j) = norm(obj_coords(id1,:)-obj_coords(neighbors_of_same_type(j),:));
               %adj_dist(id1, adj_obj(i,6+j)) = norm(obj_coords(id1,:),obj_coords(adj_obj(i,6+j),:));
               %adj_dist(adj_obj(i,6+j),id1) = adj_dist(id1, adj_obj(i,6+j));
           end
           % limit the number of neighbor to be fewer than max_num_neighbors
           num_neighbors_of_same_type = min(max_num_neighbors, num_neighbors_of_same_type);       
           [sort_dists, sort_indx] = sort(distances);
           indx_1{i} = repmat(id1-first_indx+1,[1, num_neighbors_of_same_type]);
           indx_2{i} = neighbors_of_same_type(sort_indx(1:num_neighbors_of_same_type))-first_indx+1;
           dist_values{i} = sort_dists(1:num_neighbors_of_same_type); %repmat(distances,[1, 2]);
       end
       
   end
end
%T = toc; fprintf('Indexing done in %.2f seconds\n',T);
% Elapsed time is 0.411493 seconds.
indx_1 = cat(2,indx_1{:});
indx_2 = cat(2,indx_2{:});
dist_values = cat(2,dist_values{:});
%proportion_features = cat(1, proportion_features{:});
% % sample neighboring pixels for spatial relationships
% nSamples = 10000; 
% ind_sample = randperm(length(indx_1),nSamples);
% neighbor_1 = indx_1(ind_sample);
% neighbor_2 = indx_2(ind_sample);
% figure; ndhist(proportion_features(neighbor_1,1),proportion_features(neighbor_2,1),'axis',[0 1 0 1],'filter','bins',1,'columns');
% figure; ndhist(proportion_features(neighbor_1,2),proportion_features(neighbor_2,2),'axis',[0 1 0 1],'filter','bins',1,'columns');
%adj_matrix = sparse([indx_1 indx_2],[indx_2 indx_1],1,num_total_objects,num_total_objects);
%adj_dist = sparse([indx_1,indx_2],[indx_2 indx_1],[dist_values dist_values],num_total_objects,num_total_objects);
med_dist = median(unique(dist_values)); 
% save the distance values
%save(fullfile(outdir,[imname '_distances_to_nn.mat']),'dist_values');
if obj_type == 3
    indx = ones(length(dist_values)); %dist_values >= 3*med_dist;%prctile(dist_values,75);
else
    indx = dist_values >= med_dist;%prctile(dist_values,75);
end
sigma = 1.5*med_dist;
dist_new = exp(- dist_values.^2./(2*sigma^2));
dist_new(indx) = 0;
%similarities = sparse(indx_1,indx_2,dist_new,num_total_objects,num_total_objects);
%similarities = sparse([indx_1,indx_2],[indx_2 indx_1],[dist_new dist_new],num_total_objects,num_total_objects);
similarities = sparse([indx_1,indx_2],[indx_2 indx_1],[dist_new dist_new],num_superpx,num_superpx);
%tic;
[~,C] = graphconncomp(similarities);
%T = toc; fprintf('Connected component done in %.2f seconds\n',T);

components_area = zeros(max(C),1);
for i = 1:length(components_area)
    indx_cl = C == i;
    components_area(i) = sum(obj_radii(indx_cl).^2);
end

[~,indx] = sort(components_area,'descend');
top_centers = cell(num_comps,1);
top_radii = cell(num_comps,1);
for i = 1:num_comps
   indx_cl = C == indx(i);
   top_centers{i} = adj_obj_coords(indx_cl,:);
   top_radii{i} = adj_obj_radii(indx_cl);
end
% 
% h = histogram(C,unique(C));
% num_elts = h.Values;close all; %grpstats(C,C,{'numel'});
% [~, indx] = sort(num_elts,'descend');
% num_comps = max(num_comps,sum(num_elts > 35));
% top_centers = cell(num_comps,1);
% top_radii = cell(num_comps,1);
% for i = 1:num_comps
%    indx_cl = C == indx(i);
%    top_centers{i} = obj_coords(indx_cl,:);
%    top_radii{i} = obj_radii(indx_cl);
% end
%output_dir = 'Z:\ADH_Jeff\connected_component_sparse_Feb25';
if plot_flag
    I = imread(fullfile(IMG_DIR,'images',[imname '.tif']));
    seg = zeros(nrow,ncol);
    for i = 1:num_comps
        x = top_centers{i}(:,1); y = top_centers{i}(:,2);
        if length(x) == 1
            ang=0:0.01:2*pi;
            xp=floor(top_radii{i}*cos(ang)'+x);
            yp=floor(top_radii{i}*sin(ang)'+y);
            k = boundary(xp,yp);
            mask = poly2mask(xp(k),yp(k),nrow, ncol);
        else
            k = boundary(x,y,0.5);
            mask = poly2mask(x(k),y(k),nrow, ncol);
        end
        seg = seg + mask*i;

    end
    %I2 = cat(3,redc, greenc, bluec);
    %imwrite(I2,fullfile(outdir,[imname '_obj_type_' num2str(obj_type) '.tif']));
    seg = seg + 1;
    segs{1} = seg(1:4:end, 1:4:end); % downsize the input for eval
    parsave(fullfile(outdir,[imname '.mat']),segs);

end

end

%plot(x(k),y(k),'-','Color',colors(i,:),'LineWidth',3);
%fig = figure; imshow(uint8(I)); hold on
%plot(x(k),y(k),'-k','LineWidth',3);
%hold off;
%resizeImageFig(fig, size(I(:,:,1)),2);
%saveas(fig,[outfilename '_' num2str(i) '.jpg']);
%close all;
%         mask_obj = zeros(nrow,ncol);
%         indx_obj = sub2ind(size(mask_obj),y,x);
%         mask_obj(indx_obj) = 1;
%         mask_obj = imdilate(mask_obj,se)>0;
%         %areas(i) = (max(y) - min(y))*(max(x) - min(x));
%         area_mask = sum(mask(:));
%         if area_mask > 20000
%             bdry = edge(mask);
%             bdry = imdilate(bdry,se);
%             redc(bdry) = colors(i,1);
%             redc(mask_obj) = colors(i,1);
%             greenc(bdry) = colors(i,2);
%             greenc(mask_obj) = colors(i,2);
%             bluec(bdry) = colors(i,3);
%             bluec(mask_obj) = colors(i,3);
%             %viscircles(top_centers{i}, top_radii{i},'EdgeColor',colors(i,:));
%             %plot(x(k),y(k),'-','Color',colors(i,:),'LineWidth',3);
%             %save([outfilename '_' num2str(i) '.mat'],'x','y','k','area_mask');
%         end
% %% euclidian distance between centers
% centers = obj_coords(indx_type,:);
% radii = obj_radii(indx_type);
%
% D = pdist(centers,'euclidean');
% squareformD = squareform(D);
% % sort distances
% [sortedDist, indx_dist] = sort(squareformD,1,'ascend');
% % the top row is all 0, self-self distances
% nearest_neighbor_distance = sortedDist(2:1+num_neighbors,:);
% med_dist = median(nearest_neighbor_distance(:));
% % distance threshold to identify connections between objects/superpixels
% dist_thres = med_dist; % pixels for now
% indx = D > dist_thres;
% D( indx ) = 0;
% sigma = 1.5*med_dist; % standard deviation of the distance or med*1.5
% similarities = exp(- D.^2./(2*sigma^2));
% similarities(indx) = 0;
% A = squareform(similarities);
%
% %identify connected component
% [S,C] = graphconncomp(sparse(A));
% % sort components based on number of elements
% num_elts = grpstats(C,C,{'numel'});
% [~, indx] = sort(num_elts,'descend');
% top_centers = cell(num_comps,1);
% top_radii = cell(num_comps,1);
% for i = 1:num_comps
%    indx_cl = C == indx(i);
%    top_centers{i} = centers(indx_cl,:);
%    top_radii{i} = radii(indx_cl);
% end

%n_n_list = unique(n_n(:)); n_n_list(n_n_list == 0 | n_n_list == id1) = [];
%neighbors of neighbors of neighbors - MUTUAL KNN
%n_n_n = adj_map(n_n_list',[1,7:end]);
% find mutual knn (both are neighbors of each other)
%for j = 1:size(n_n_n,1)
%    if sum(n_n_n(j,:) == id1) == 0
%        n_n_n(j,:) = 0;
%    end
%end
%neighbor_indx = n_n_n(:,1)';neighbor_indx(neighbor_indx == 0) = [];
       
       %neighbors of neighbors of neighbors
       %n_n_n = adj_map(n_n,7:end);n_n_n = n_n_n(:)';
       %n_n_n(n_n_n == 0) = [];
       %neighbor_indx = [adj_obj(i,7:(6+num_neighbors_input)), n_n, n_n_n];       
       
          %I = double(I)./255; redc = I(:,:,1); greenc = I(:,:,2); bluec = I(:,:,3);
    %outfilename = fullfile(outdir,[imname param_string '_conncomp_obj_type' num2str(obj_type)]);
    %colors = distinguishable_colors(num_comps,[0 1 0]);
    %se = strel('disk',7,4);


