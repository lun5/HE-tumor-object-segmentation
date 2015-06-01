
labels = bwlabel(E_ucm < thresh);
numRegions = length(unique(labels))-1;
% we have affinity_matrix
% load local pair
cache_file = sprintf('caches/ii_jj_caches/%d_%d_%d.mat',im_sizes{1}(1),im_sizes{1}(2),7);
data = load(cache_file); data = data.data;
ii = data.ii;
jj = data.jj;
W = affinity_matrix{1}; thr_affinity = prctile(nonzeros(W),75);
Area = im_sizes{1}(1)*im_sizes{1}(2);
for i = 1:numRegions - 1
    for j = i+1:numRegions
        reg1 = labels == i;
        reg2 = labels == j;
        if sum(reg1(:)) == 0 || sum(reg2(:)) == 0 || sum(reg1(:)) > Area/16 || sum(reg2(:)) > Area/16
            continue;
        end
        [reg1_row, reg1_col] = find(reg1);
        [reg2_row, reg2_col] = find(reg2);
        if max(reg1_row) < min(reg2_row) || min(reg1_row) > max(reg2_row) ||...
                max(reg1_col) < min(reg2_col) || min(reg1_col) > max(reg2_col) 
            continue; % totally out of each other range
        end
        ind_reg1 = find(reg1);
        ind_reg2 = find(reg2);
        ind_reg1_in_ii = ismember(ii,ind_reg1);
        ind_reg1 = ii(ind_reg1_in_ii);
        ind_neigh_reg1 = jj(ind_reg1_in_ii);
        ind_neigh_reg1_in_reg2 = ismember(ind_neigh_reg1,ind_reg2);
        if sum(ind_neigh_reg1_in_reg2) == 0
            continue; 
        end;
        ind_1 = ind_reg1(ind_neigh_reg1_in_reg2);
        ind_2 = ind_neigh_reg1(ind_neigh_reg1_in_reg2);
        ind_neighboring_pixels = sub2ind(size(W),ind_1,ind_2);
        affinity_reg1_reg2 = W(ind_neighboring_pixels);
        if sum(affinity_reg1_reg2 > thr_affinity) / size(affinity_reg1_reg2,1) > 1/4
            fprintf('Merge Region %d and %d\n',i,j);
            labels(labels == j) = i;
        end
    end
end
pause;