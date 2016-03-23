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
%tiles_dir = fullfile('Z:\','TilesForLabeling_tiff_renamed'); %window
%tiles_dir = '/Users/lun5/Research/data/TilesForLabeling_tiff_renamed'; %mac
%tiles_dir = '/home/lun5/HEproject/data/normalization_512';
%tiles_dir = '/Users/lun5/Research/data/normalization_2048_tif';
%tiles_dir = '/Users/lun5/Research/data/normalization_512';
%tiles_dir = '/Users/lun5/Research/data/normalization_512';
%tiles_dir =  '/home/lun5/HEproject/TilesForLabeling_tiff_renamed'; %linux
%tiles_dir = 'Z:\HEproject\normalization_512/invasive';
tiles_dir = 'Z:\HEproject\normalization_512/well_defined';
clear raw_image Pts ans im mdist opts_affinity opts_clustering which_affinity which_features
%imname = 'djuten6dhnfd.tif';
%imname = '4d0ylpdlwf.tif';
imname = '8ghygsmwjy.tif';
%imname = 'hrlxfimwgjas.tif';
%imname = 'uaZFwoHref.tif';
%imname = 'aNaggwovpxANWq0.tif';
%imname = 'jRh62FQ8hUZWlA.tif';
%imname = '0ANZqyIBfUc.tif';
%imname ='aW5aZV9o5NgqVX.tif';
%imname = '95f7k8loesyevi.tif';
%imname = 'cxwrYBYWredN.tif';
%imname = 'JDXGoRjONolk.tif';
%imname = 'p76gode3evdqoin.tif';
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
%imname = 'dj0ebjbuxshxz.tif';
%imname = 'aqizfuqbbxyu.tif';
%imname = 'uraxeh1spli7ky9.tif';
imname = 'q9VDQzxnxb.tif';
%imname = '1BHJQxeCcT.tif';
%imname = 'TaLJYO23jlXd.tif';
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
%imname = 'P76gODe3EVdqOiN.tif';
%imname ='klmqi6sq7wl6.tif'
%% result directory
% splitStr = regexp(imname,'\.','split');
% imresult_dir = fullfile(pwd,'results','HE_results',[splitStr{1} 'crop2']);
% 
% if ~exist(imresult_dir,'dir')
%     mkdir(imresult_dir);
%     fileattrib(imresult_dir,'+w');
% end
imname = lower(imname);
%imname = fullfile(pwd,'results','im.tif');dz_im = imread(imname);
raw_image = imread(fullfile(tiles_dir, imname));%figure;imshow(raw_image);
ndown = 1;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
%figure; imshow(raw_image);rect = getrect;dz_im = imcrop(raw_image,rect);
I = double(dz_im);figure; imshow(I/255);size(I)
%I = double(raw_image);
%%
% Pts array is updated
%I = raw_image;
opts_affinity = setEnvironment_affinity;
opts_affinity.features.which_features = {'hue opp'};%, 'brightness opp', 'saturation opp'};
opts_affinity.joint_exponent = 1.5; opts_affinity.sig = 3; %opts.p_reg = 100;
[Ws,Ws_each_feature_set, im_sizes] = getW(I,opts_affinity);
% %% Graph-based clustering based on
% % this depends on whether the outputs are segmentation or detecting edges
opts_clustering = setEnvironment_clustering;
[segmented_image_allfeatures,E_ucm_weighted, E_weighted, E_oriented] = graphSegmentation(Ws,Ws_each_feature_set{1},im_sizes,I,opts_clustering);
disp('Done');
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