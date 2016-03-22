% script for testing SIC
% Luong Nguyen 3/20/2016

clearvars; close all;
tiles_dir = '/Users/lun5/Research/data/Tiles_512';
%imname = 'k2yxq1tbr6kpny0.tif';imname = 'aNaggwovpxANWq0.tif';imname ='klmqi6sq7wl6.tif';
%imname = 'P76gODe3EVdqOiN.tif';imname = 'TaLJYO23jlXd.tif';imname = '7xn9dazygam.tif';
%imname = 'dRfMkOErZY.tif';%imname = 'fFwTGXYlhYNa.tif';imname = 'ibhyyugefpbn.tif';
%imname = 'lszomrlgsc5na4q.tif';imname = 'uraxeh1spli7ky9.tif';imname = 'vm3qo9caekfodi.tif';
%imname = '95f7k8loesyevi.tif';imname = 'JDXGoRjONolk.tif';
%imname = 'w8kwtop6hyp.tif';imname = 'h1402uhfkz.tif';imname = 'ycivjoxn14stvq.tif';
%imname = 'mbdqhorkuxs.tif';imname = 'pLYZEV43nHWmUDK.tif';imname = 'jbaKL4TsEqT.tif';
%imname = 'h1402uhfkz.tif';imname = 'jRh62FQ8hUZWlA.tif';imname ='aW5aZV9o5NgqVX.tif';
imname = 'h1402uhfkz.tif';
imname = lower(imname);
%imname = fullfile(pwd,'results','im.tif');dz_im = imread(imname);
raw_image = imread(fullfile(tiles_dir, imname));%figure;imshow(raw_image);
ndown = 1;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
%figure; imshow(raw_image);rect = getrect;dz_im = imcrop(raw_image,rect);
im_rgb = double(dz_im);figure; imshow(im_rgb./255);

% extract random pixels to get the rotation matrix
%nrow = size(im_rgb,1); ncol = size(im_rgb,2);
% random_indx = randperm(nrow*ncol, 50000);
% X = reshape(im_rgb,[nrow*ncol 3]);
% data = X(random_indx,:);
% mean_im = mean(data,1);
% mean_im = mean_im./norm(mean_im);
% mean_comp = data*mean_im';
% left_over = data - mean_comp*mean_im;
% [u, s, v] = svd(left_over,'econ');
% rotation_mat = v';
% rotation_mat = cat(1,abs(rotation_mat(3,:)), rotation_mat(1:2,:));
rgb_coords = reshape(im_rgb,[size(im_rgb,1)*size(im_rgb,2),size(im_rgb,3)])';
rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
rotation_matrix = rotation_matrix.rotation_matrix;
mu_s = .1; sigma_s = 1;
[ sic_coords ] = rgb2sic( rgb_coords, mu_s, sigma_s, rotation_matrix);
sic1_im = reshape(sic_coords(1,:),[size(im_rgb,1), size(im_rgb,2)]);
sic2_im = reshape(sic_coords(2,:),[size(im_rgb,1), size(im_rgb,2)]);
figure; imagesc(sic1_im); axis off; axis square; colorbar; caxis([-1 1]);
figure; imagesc(sic2_im); axis off; axis square; colorbar; caxis([-1 1]);

nrows = size(im_rgb,1); ncols = size(im_rgb,2);
num_pixels = nrows*ncols;
Nsamples = 10000;
opts = setEnvironment_affinity;
opts.sig = 3;
f_maps_curr = sic1_im;
F = sampleF(f_maps_curr,Nsamples,opts);

figure; ndhist(F(:,1),F(:,2),'axis',[-1 1.2 -1 1.2],'filter','bins',1,'columns');
xlabel('sic1_A'); ylabel('sic1_B');
axis square;set(gca,'FontSize',16);set(gcf,'color','white')

% X = f_maps_curr(:);
% GMModel1 = fitgmdist(X,2,'Start','plus');
% 
% Mu = [GMModel1.mu(1) GMModel1.mu(1); GMModel1.mu(1) GMModel1.mu(2);...
%      GMModel1.mu(2) GMModel1.mu(1);GMModel1.mu(2) GMModel1.mu(2)];
% Sigma = zeros(2,2,4);
% Sigma(:,:,1) = diag([1 1]);
% Sigma(:,:,2) = diag([1 1]);
% Sigma(:,:,3) = diag([1 1]);
% Sigma(:,:,4) = diag([1 1]);
% 
% PComponents = [1/4,1/4,1/4,1/4];
% S = struct('mu',Mu, 'Sigma',Sigma,'ComponentProportion',PComponents);
% GMModel3 = fitgmdist(F,4,'Start',S);
% 
% mu = GMModel3.mu;
% sigma = GMModel3.Sigma;
% prop = GMModel3.ComponentProportion;
% obj = gmdistribution(mu,sigma,prop);

% ezsurf(@(x,y)pdf(obj,[x y]),[-1 1.2],[-1 1.2])

%% using KDE
opts.features.which_features = {'hue sat'};
which_features = {'hue sat'};
f_maps(:,:,1) = sic1_im;
%f_maps(:,:,2) = sic2_im;
%f_maps(:,:,3) = mean(im_rgb,3)./255;
%f_maps = f_maps(1:2:end,1:2:end,:);
p{1} = learnP_A_B(f_maps,opts);
%figure; imagesc(hist(p{1})); set(gca,'FontSize',20); axis ij; axis xy;axis square; colorbar
rf{1} = learnPMIPredictor(f_maps,p{1},[], 'hue-sat', opts);

Ws = buildW_pmi(f_maps,rf{1},p{1},[], which_features, opts);

opts_clustering = setEnvironment_clustering;
[~,~,~, E_oriented] = graphSegmentation([],{Ws},{size(f_maps(:,:,1))},im_rgb,opts_clustering);
