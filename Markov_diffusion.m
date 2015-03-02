sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%tiles_dir = '/Users/lun5/Box Sync/TilesForLabeling_tiff_renamed';
opts_input = setEnvironment_inputs;

I = imread(fullfile(tiles_dir, '9uixINHtjjiS.tif'));
I = imresize(I,1/4);
imshow(I);
rect = getrect; rect = round(rect);
I = imcrop(I,rect);imshow(I);size(I)
%I = imread(fullfile(pwd,'test_images','253027.jpg'));
%imwrite(I,fullfile('test_images','tp10-611gland7snip.tif'),'tif','Compression','none');
%% Calculate affinity matrix 
%I_downsample = imresize(I,1/4);figure;imshow(I_downsample);
opts_affinity = setEnvironment_affinity;
[affinity_matrix, im_sizes] = calculateAffinity(I, opts_affinity);

orig_sz = size(I);
tx = orig_sz(1);
ty = orig_sz(2);

W = affinity_matrix{1};
[wx, wy] = size(W);
x = 1 : wx;
S = full(sum(W, 1));
D = sparse(x, x, S, wx, wy);
M = W/D; % markov transition matrix

opts.issym=1;
opts.isreal = 1;
opts.disp=0;
nvecs = 20;
[Evecs, Evals] = eigs(W/D.^(1/2), D.^(1/2), nvecs, 'sm',opts); %
U = Evecs;
dm = 6;
foptions(1) = 1;foptions(14) = 0;
foptions(2) = 1e-2; foptions(3) = 1e-2;
tstKmeansSpectral;

[~,indx_membership] = max( post > 0,[],2);
for cl = 1:6
    id_cluster = reshape(indx_membership, [orig_sz(1) orig_sz(2)]);
    id_cluster(id_cluster ~=cl) = 0;
    id_cluster(id_cluster ~=0) = 1;
    id_im = uint8(I).*uint8(repmat(id_cluster,1,1,3));
    h=figure; imshow(id_im);
    set(gcf,'color','white') % White background for the figure.
end
start_pos = 589; num_iter = 5;
%state_vecs = zeros(orig_sz(1)*orig_sz(2),num_iter);
%state_vecs(start_pos,1) = 1;
transition_matrix = M;%speye(orig_sz(1)*orig_sz(2)); 
%F(num_iter -1) = struct('cdata',[],'colormap',gray);
state_im = zeros(orig_sz(1),orig_sz(2),num_iter-1);
for iter = 1:3%num_iter-1
   transition_matrix = transition_matrix*transition_matrix;
   prob_vec = transition_matrix(:,start_pos); % probability at next time step
   state_im(:,:,iter) = reshape(prob_vec,[orig_sz(1) orig_sz(2)]);
   [X,map] = gray2ind(state_im(:,:,iter));
   %figure; imshow(X);
   imtool(log(state_im(:,:,iter)+eps),[])
   %F(iter) = im2frame(X,map);
end
% %movie2avi(F,'MarkovMovie.avi','compression','None');
% %movie(F,2);
% %% Graph-based clustering based on 
% % this depends on whether the outputs are segmentation or detecting edges
opts_clustering = setEnvironment_clustering;
[segmented_image, E_oriented] = graphSegmentation(affinity_matrix,im_sizes,I,opts_clustering,opts_affinity);
