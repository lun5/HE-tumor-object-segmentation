%sourcedir = 'Z:\';
%tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%tiles_dir = '/Users/lun5/Box Sync/TilesForLabeling_tiff_renamed';
tiles_dir = fullfile(pwd,'HEimages');
%opts_input = setEnvironment_inputs;

raw_image = imread(fullfile(tiles_dir, '9uixINHtjjiS.tif'));
%raw_image = imread(fullfile(tiles_dir, 'EMnOxgxqoMGzn1.tif'));
%raw_image = imresize(raw_image,1/4);
figure;imshow(raw_image);
rect = getrect; rect = round(rect);
raw_image = imcrop(raw_image,rect);imshow(raw_image);size(raw_image)
I = double(raw_image);
%I = imread(fullfile(pwd,'test_images','random206863.pgm'));
%%
% Pts array is updated
opts_affinity = setEnvironment_affinity;
tic;
[Pts,A,mdist] = calculateAffinity(I, opts_affinity);
disp('fast calculation?');toc
% A = affinity_matrix{1};

sizeIm = size(I(:,:,1));
%Pts = ones(prod(sizeIm),2);
%Pts(:,1) = I(:);
% afftyPar.sizeIm  = sizeIm;
% afftyPar.dsThres = 3.1; %kanisza
% afftyPar.rho     = 1.5; %kanisza
% tic;
% [Pts,A1,binNhbr,mdist1] = mkAffty(Pts,afftyPar);
% disp('slow calculation?');toc
[discComp,U,S] = clusterKmeans(A,4);
displayComp(Pts,discComp,1,sizeIm)

% [discComp,U,S] = clusterKmeans(A1,4);
% displayComp(Pts*255,discComp,1,sizeIm)

%imwrite(I,fullfile('test_images','tp10-611gland7snip.tif'),'tif','Compression','none');
%% Calculate affinity matrix 

vect = zeros(sizeIm(1),sizeIm(2),size(U,2));
for i =1:size(U,2);
    vect(:,:,i) = reshape(U(:,i),sizeIm(1),sizeIm(2));
end

figure;
ha = tight_subplot(2,2,[.01 .0],[0 0],[0 0]);
for i =1:4
    axes(ha(i));imagesc(vect(:,:,i));
    axis equal; axis tight; axis off; %colormap('gray');
end
set(gcf,'color','white') 


% D = sum(A, 1)';              % Normalize column sum to one.
% D = sparse(1:prod(sizeIm),1:prod(sizeIm),D,prod(sizeIm),prod(sizeIm));
% M = A/D; % markov transition matrix
% start_pos = 4; num_iter = 5;
% transition_matrix = M;
% state_im = zeros(sizeIm(1), sizeIm(2),num_iter-1);
% for iter = 1:num_iter-1
%    transition_matrix = transition_matrix*transition_matrix;
%    prob_vec = transition_matrix(:,start_pos); % probability at next time step
%    state_im(:,:,iter) = reshape(prob_vec,[sizeIm(1) sizeIm(2)]);
%    imtool(state_im(:,:,iter),[])
% end

% %% Graph-based clustering based on
% % this depends on whether the outputs are segmentation or detecting edges
opts_clustering = setEnvironment_clustering;
[segmented_image, E_oriented] = graphSegmentation({A},{[sizeIm(2) sizeIm(1)]},I,opts_clustering,opts_affinity);
