%sourcedir = 'Z:\';
%tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%tiles_dir = '/Users/lun5/Box Sync/TilesForLabeling_tiff_renamed';
tiles_dir = fullfile(pwd,'HEimages');

%imname = '9uixINHtjjiS.tif';
imname = '2ALe5NgRyfnpo.tif';
%% result directory
splitStr = regexp(imname,'\.','split');
imresult_dir = fullfile(pwd,'results','HE_results',splitStr{1});

if ~exist(imresult_dir,'dir')
    mkdir(imresult_dir);
    fileattrib(imresult_dir,'+w');
end

% raw_image = imread(fullfile(tiles_dir, imname));
% %raw_image = imread(fullfile(tiles_dir, 'EMnOxgxqoMGzn1.tif'));
% % raw_image = imresize_local(raw_image,3);
% image(raw_image); axis off; axis equal;
% rect = getrect;%[919.551244509517 580.716691068814 152.92532942899 113.944363103953];%getrect; 
% rect = round(rect);
% crop_image = imcrop(raw_image,rect);size(crop_image)
% figure;image(crop_image); axis off; axis equal;
% I = double(crop_image);
% 
% save the original image
imwrite(crop_image,fullfile(imresult_dir,'crop_image.tif'));
% I = double(crop_image.*255);
%I = double(raw_image*255);
%I = imread(fullfile(pwd,'test_images','random206863.pgm'));
%%
% Pts array is updated
%I = raw_image;
opts_affinity = setEnvironment_affinity;
which_features = opts_affinity.features.which_features;
which_affinity = opts_affinity.affinityFunction;
methodresult_dir = fullfile(imresult_dir,[which_features{1} '_' which_affinity]);
if ~exist(methodresult_dir,'dir')
    mkdir(methodresult_dir);
    fileattrib(methodresult_dir,'+w');
end

tic;
[Pts,A,mdist] = calculateAffinity(I, opts_affinity);
disp('fast calculation?');toc
%imStats(full(A));
% A = affinity_matrix{1};
sizeIm = size(I(:,:,1));
im = reshape(Pts,sizeIm);
% %% Graph-based clustering based on
% % this depends on whether the outputs are segmentation or detecting edges
opts_clustering = setEnvironment_clustering;
[segmented_image, E_oriented] = graphSegmentation({A},{[sizeIm(2) sizeIm(1)]},I,opts_clustering,opts_affinity);
parsave(fullfile(methodresult_dir,'E_oriented'),E_oriented);


%Pts = ones(prod(sizeIm),2);
%Pts(:,1) = I(:);
% afftyPar.sizeIm  = sizeIm;
% afftyPar.dsThres = 3.1; %kanisza
% afftyPar.rho     = 1.5; %kanisza
% tic;
% [Pts,A1,binNhbr,mdist1] = mkAffty(Pts,afftyPar);
% disp('slow calculation?');toc
% tic;[discComp,U,S] = clusterKmeans(A,4);toc
% displayComp(Pts,discComp,1,sizeIm(1:2))
% 
% %% Calculate affinity matrix 
% 
% vect = zeros(sizeIm(1),sizeIm(2),size(U,2));
% for i =1:size(U,2);
%     vect(:,:,i) = reshape(U(:,i),sizeIm(1),sizeIm(2));
% end
% 
% figure;
% ha = tight_subplot(2,2,[.01 .0],[0 0],[0 0]);
% for i =1:4
%     axes(ha(i));imagesc(vect(:,:,i));
%     axis equal; axis tight; axis off; %colormap('gray');
% end
% set(gcf,'color','white') 

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

