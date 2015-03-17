aff_matrix = Ws{1};
im_theta = f_maps{1};
indx = find(aff_matrix(1,:))
im_theta(1:4,1:4)
reshape(full(aff_matrix(1,indx)),[4 4])

aff_matrix = affinity_matrix{1};
im_theta = f_maps{1};
indx = find(aff_matrix(516,:))
im_theta(1:4,1:4)
reshape(full(aff_matrix(1,indx)),[4 4])
