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
%imname = '4d0ylpdlwf.tif';
%imname = '8ghygsmwjy.tif';
%imname = 'hrlxfimwgjas.tif';
%imname = 'uaZFwoHref.tif';
%imname = 'aNaggwovpxANWq0.tif';
%imname = 'jRh62FQ8hUZWlA.tif';
%imname = '0ANZqyIBfUc.tif';
imname = '95f7k8loesyevi.tif';
%imname = 'cxwrYBYWredN.tif';
%imname = '9uixINHtjjiS.tif';
%imname = 'w8kwtop6hyp.tif';
%imname = '2ALe5NgRyfnpo.tif';
%imname = 'jbaKL4TsEqT.tif';
%imname = 'k2yxq1tbr6kpny0.tif';
%imname = 'vmp8mdxkod3xxzu.tif';
%imname = 'vm3qo9caekfodi.tif';
%imname = 'h1402uhfkz.tif';
%imname = 'dRfMkOErZY.tif';
%imname = 'ycivjoxn14stvq.tif';
%imname = 'fFwTGXYlhYNa.tif';
%imname = 'ibhyyugefpbn.tif';
%imname = 'pLYZEV43nHWmUDK.tif';
%imname = '7xn9dazygam.tif';
%imname = 'mbdqhorkuxs.tif';
%imname = 'lszomrlgsc5na4q.tif';
%imname = 'n2wolhpsak70anw.tif';
%imname = '7vj4ekusieek6ys.tif';
%imname = 'aqizfuqbbxyu.tif';
% imnames = {'aux48hgyn767ebt.tif','ervrkyrmpb.tif','finkidqlnznihk.tif',...
%     'fs4dqvvb2xu.tif','fs5tuzogn0.tif','gg9wmyahfpc0c.tif','hrdqmlu2ig.tif',...
%     'iitv3ixlbhih5q.tif','k3gutsibfg1qxg7.tif','klmqi6sq7wl6.tif','lc20hzhj6p.tif',...
%     'lszomrlgsc5na4q.tif','n2wolhpsak70anw.tif','oasyczs5983.tif','p1edtnwreykh.tif',...
%     'rnky9htjaz.tif','s5rkUDO7teXwKl.tif','sfZyWlg4291cH9.tif','shyY81fuot9HU.tif',...
%     'sp2BYg1b33Ghl.tif','TaLJYO23jlXd.tif','UHIN9NL4Ju7BLS.tif','UYUuKDfZQq.tif',...
%     'wM9G9oyReU.tif','XHCY7eRoyqZ.tif','xnnGATe0Qi1mL.tif','YLf5dGhxYp2Sxp.tif',...
%     'ZcU7vNvnWSrKAf.tif','zhB3vhJ2vc.tif','ZleLRgfLYYfT4C.tif','zpD2rLbjZm.tif','zXvCwQJOEyD.tif'};
%i = 15;
%imname = '13nedzdzfh.tif';
%imname = imnames{i}
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
ndown = 4;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
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
%[A,im_sizes] = getW(I,opts_affinity);
[Ws,Ws_each_feature_set, im_sizes] = getW(I,opts_affinity);
%disp('fast calculation?');toc

%sizeIm = size(I(:,:,1));
%im = reshape(Pts,sizeIm);
% %% Graph-based clustering based on
% % this depends on whether the outputs are segmentation or detecting edges
opts_clustering = setEnvironment_clustering;
[segmented_image, E, E_oriented] = graphSegmentation(Ws,Ws_each_feature_set{1},im_sizes,I,opts_clustering);
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