aff_matrix = Ws{1};
im_theta = f_maps{1};
indx = find(aff_matrix(1,:))
im_theta(1:4,1:4)
reshape(full(aff_matrix(1,indx)),[4 4])

    -0.3785   -0.4567   -0.3697   -2.4201
   -1.4555   -0.8981    0.0578   -0.3286
   -0.2119   -0.0380    2.5055    0.1947
   -0.1117    0.0017    0.5766   -0.3662
   
    1.0000    0.9986    1.0000    0.1389
    0.1422    0.6489    0.9952    0.9997
    0.9979    0.9958    0.1355    0.9946
    0.9965    0.9955    0.9896    1.0000

aff_matrix = affinity_matrix{1};
im_theta = f_maps{1};
indx = find(aff_matrix(516,:))
im_theta(1:4,1:4)
reshape(full(aff_matrix(1,indx)),[4 4])
