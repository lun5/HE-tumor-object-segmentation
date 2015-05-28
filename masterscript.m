% Demo file
% -------------------------------------------------------------------------
% HE tumor object segmentation 
% Luong Nguyen, lun5@pitt.edu
% Please email me if you find bugs, or have suggestions or questions

function masterscript

%% compile and check for error
%addpath(genpath(pwd));
%compile;

%%
%tiles_dir = fullfile(pwd,'data','images','test');
%sourcedir = 'Z:\';
tiles_dir = fullfile('Z:\','TilesForLabeling_tiff_renamed'); %window
%tiles_dir = '/Users/lun5/Box Sync/TilesForLabeling_tiff_renamed'; %mac
%tiles_dir =  '/home/lun5/HEproject/TilesForLabeling_tiff_renamed'; %linux
clear raw_image Pts ans im mdist opts_affinity opts_clustering which_affinity which_features
%imname = 'uaZFwoHref.tif';
%imname = 'aNaggwovpxANWq0.tif';
%imname = 'jRh62FQ8hUZWlA.tif';
%imname = '0ANZqyIBfUc.tif';
%imname = '95f7k8loesyevi.tif';
%imname = 'cxwrYBYWredN.tif';
%imname = '9uixINHtjjiS.tif';
%imname = 'w8kwtop6hyp.tif';
%imname = '2ALe5NgRyfnpo.tif';
%imname = 'jbaKL4TsEqT.tif';
%imname = 'k2yxq1tbr6kpny0.tif';
%imname = 'vmp8mdxkod3xxzu.tif';
%imname = 'h1402uhfkz.tif';
imname = 'dRfMkOErZY.tif';
%imname = 'ycivjoxn14stvq.tif';
%imname = 'fFwTGXYlhYNa.tif';
%imname = 'pLYZEV43nHWmUDK.tif';
%imname = 'LLV232_D04_20x_max_proj.tif';
%tiles_dir = fullfile(pwd,'test_images');
%imname = '253027.jpg';
%% result directory
% splitStr = regexp(imname,'\.','split');
% imresult_dir = fullfile(pwd,'results','HE_results',[splitStr{1} 'crop2']);
% 
% if ~exist(imresult_dir,'dir')
%     mkdir(imresult_dir);
%     fileattrib(imresult_dir,'+w');
% end
imname = lower(imname);
raw_image = imread(fullfile(tiles_dir, imname));
ndown = 2;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
%figure; imshow(raw_image);rect = getrect;dz_im = imcrop(raw_image,rect);
I = double(dz_im);figure; imshow(I/255);size(I)
%I = double(raw_image);
%%
% Pts array is updated
%I = raw_image;
opts_affinity = setEnvironment_affinity;
%which_features = opts_affinity.features.which_features;
%which_affinity = opts_affinity.affinityFunction;
% methodresult_dir = fullfile(imresult_dir,[which_features{1} '_' which_affinity]);
% if ~exist(methodresult_dir,'dir')
%     mkdir(methodresult_dir);
%     fileattrib(methodresult_dir,'+w');
% end

%tic;
%[Pts,A,mdist] = calculateAffinity(I, opts_affinity);
[A,im_sizes] = getW(I,opts_affinity);
%disp('fast calculation?');toc

%sizeIm = size(I(:,:,1));
%im = reshape(Pts,sizeIm);
% %% Graph-based clustering based on
% % this depends on whether the outputs are segmentation or detecting edges
opts_clustering = setEnvironment_clustering;
[segmented_image, E, E_oriented] = graphSegmentation(A,im_sizes,I,opts_clustering);
%parsave(fullfile(methodresult_dir,'E_oriented'),E_oriented);

end
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