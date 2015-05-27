% how to show that the joint distribution is unique
% Given an image with some structure
% randomize the image
% look for the changes in the joint distribution
% this could be shown by fitting the joint distribution
% or this could be shown by the raw 

% load image
close all; clearvars;
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%imname = '9uixINHtjjiS.tif';
%imname = '2ALe5NgRyfnpo.tif';
%imname = 'jbaKL4TsEqT.tif';
%imname = 'k2yxq1TBR6kpNY0.tif';
%imname = 'jRh62FQ8hUZWlA.tif';
imname = 'dRfMkOErZY.tif';
%imname = 'fFwTGXYlhYNa.tif';
%imname = 'pLYZEV43nHWmUDK.tif';

splitStr = regexp(imname,'\.','split');
raw_image = imread(fullfile(tiles_dir, imname));
I = double(raw_image); 
figure; imshow(raw_image);axis tight;
% calculate the feature
opts_affinity = setEnvironment_affinity;
which_features = opts_affinity.features.which_features;
f_maps = getFeatures(double(I),1,which_features,opts_affinity);
f_map_opphue = f_maps{1};
% sample from the image
Nsamples = 10000;
F = sampleF(f_map_opphue,Nsamples,opts_affinity);  
figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
zlim([0 0.05])
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 

% randomize pixels in the image
ind_row = randperm(size(I,1));
ind_col = randperm(size(I,2));
randomized_image = raw_image(ind_row, ind_col,:);
figure; imshow(randomized_image); axis tight
% resample from the mangled image
f_map_opphue_rand = f_map_opphue(ind_row,ind_col);
% sample from the image
F = sampleF(f_map_opphue_rand,Nsamples,opts_affinity);  
figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
zlim([0 0.05])
xlabel('\phi_A'); ylabel('\phi_B'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 

%%
% fake_im = zeros(2,2,3);
% fake_im(1,1,:) = [1 1 1];
% fake_im(1,2,:) = [1 0 0];
% fake_im(2,1,:) = [0 1 0];
% fake_im(2,2,:) = [0 0 1];
% 
% figure; imshow(fake_im);
% 
% for i = 1:10
% fake_im_rand = fake_im(randperm(size(fake_im,1)),randperm(size(fake_im,2)),:);
% figure; imshow(fake_im_rand);
% end