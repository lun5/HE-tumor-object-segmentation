% Isolate white pixels by clustering the ab space into 11 clusters
% manually picked out white pixels?
% automate this: start with a group of completely white pixels (255, 255, 255)
% choose the cluster center with this images is closest to
% Luong Nguyen 3/2/2016
sourcedir = 'Z:\';
addpath(genpath(pwd));
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
%mixture_vonMises_dir = fullfile(sourcedir,'mixture_von_mises','same_rot_renamed_images_FreezeAll');

num_pix_per_tile = 20000;
num_images = 100;
fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
% rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
% rotated_coordinates = rotation_matrix.rotation_matrix*X';
% r = im_rgb(:,:,1);
% theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:)); %hue
% im_theta = reshape(theta,size(r));
% sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
% im_sat = reshape(sat,size(r));
% brightness = rotated_coordinates(1,:);
% im_brightness = reshape(brightness,size(r));
sample_ab = zeros(num_images*num_pix_per_tile,2);
sample_rgb = zeros(num_images*num_pix_per_tile,3);
cform = makecform('srgb2lab'); % transformation into Lab space
for i = 1: num_images
    tic;
    imname = imagepaths{i*2}(1:end-4)
    %imname = 'dRfMkOErZY.tif';%imagepaths{20};
    raw_image = imread(fullfile(tiles_dir, [imname '.tif']));%figure; imshow(raw_image);
    im_rgb = double(raw_image)./255;
    X = reshape(im_rgb,[size(im_rgb,1)*size(im_rgb,2),size(im_rgb,3)]);
    lab_he = applycform(im_rgb,cform);
    ab = double(lab_he(:,:,2:3));
    nrows = size(ab,1);
    ncols = size(ab,2);
    ab = reshape(ab,nrows*ncols,2);
    indx_rand = randperm(nrows*ncols,num_pix_per_tile);
    sample_ab((i-1)*num_pix_per_tile + 1:i*num_pix_per_tile,:) = ab(indx_rand,:);
    sample_rgb((i-1)*num_pix_per_tile + 1:i*num_pix_per_tile,:) = X(indx_rand,:);
    toc
end

nColors = 11;
% repeat the clustering 3 times to avoid local minima
tic;
[cluster_idx, cluster_center] = kmeans(sample_ab,nColors,'Start','plus',...
    'distance','sqEuclidean','Replicates',3);
toc   

rgb_cluster = zeros(nColors, 3);figure;
for i = 1:nColors
    indx_clusters = cluster_idx == i;
    rgb_cluster(i,:) = mean(sample_rgb(cluster_idx == i,:),1);
    rectangle('Position',[0, (i-1)*5, 5,10],'FaceColor',rgb_cluster(i,:)); hold on; 
end
axis off;
% identify which is a white cluster
white_pixel = reshape([1 1 1],[1 1 3]);
white_lab = applycform(white_pixel,cform);
ab = reshape(double(white_lab(:,:,2:3)),[1 2]);
dist_to_cluster_centers = sum((repmat(ab,[nColors 1]) - cluster_center).^2,2);
[~,white_cluster_indx] = min(dist_to_cluster_centers);
white_cluster_center = cluster_center(white_cluster_indx,:);

%% test on an image
imname = '2e95IAuSax5HD';
imname = 'dRfMkOErZY';
tic;
raw_image = imread(fullfile(tiles_dir, [imname '.tif']));%figure; imshow(raw_image);
im_rgb = double(raw_image)./255;
X = reshape(im_rgb,[size(im_rgb,1)*size(im_rgb,2),size(im_rgb,3)]);
lab_he = applycform(im_rgb,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
dist_to_clusters = zeros(nrows*ncols,nColors);
for i = 1:nColors
   dist_to_clusters(:,i) = sum((ab - repmat(cluster_center(i,:),[nrows*ncols,1])).^2,2); 
end
toc

[mm,ii] = min(dist_to_clusters,[],2);
indx_white = ii == white_cluster_indx;
XX = X;
XX(~indx_white,:) = 0;
img_out = reshape(XX,[nrows,ncols,3]); figure; imshow(img_out);

rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
rotated_coordinates = rotation_matrix.rotation_matrix*X';
r = im_rgb(:,:,1);
theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:)); %hue
mu_white = 2.24; kappa_white = 30;
theta(:,indx_white) = circ_vmrnd(mu_white, kappa_white, sum(indx_white));
im_theta = reshape(theta,size(r));
sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
sat(indx_white,:) = 0;
im_sat = reshape(sat,size(r));
brightness = rotated_coordinates(1,:);
im_brightness = reshape(brightness,size(r));

% X_cart = [rotated_coordinates(2,:); rotated_coordinates(3,:)];
X_cart = [cos(theta); sin(theta)]';
X_cart_nowhite = X_cart(~indx_white,:);
numClusters = 2;opts_mixture.noise = 1;
%% Call the function
[ mu_hat_polar,~, kappa_hat,posterior_probs, prior_probs] =...
    moVM_fixWhite(X_cart_nowhite,numClusters,opts_mixture);
%membership
[~, indx_membership] = max(posterior_probs,[],2); % 4 is the uniform noise
all_posterior_probs = zeros(ncols*nrows,4);
all_posterior_probs(~indx_white,[1:2 4]) = posterior_probs;
all_posterior_probs(indx_white,3) = 1;
post_images = reshape(all_posterior_probs,[nrows ncols 4]);

for i = 1:4
   figure; imagesc(post_images(:,:,i)); axis off; axis square
   colorbar; caxis([0 1]);
end

for cl = 1:(numClusters+1+opts_mixture.noise)
%     id_cluster = reshape(all_indx_membership, size(im_theta));
%     id_cluster(id_cluster ~=cl) = 0;
%     id_cluster(id_cluster ~=0) = 1;
%    id_im = uint8(raw_image).*uint8(repmat(id_cluster,1,1,3));
    id_im = raw_image.*uint8(repmat(post_images(:,:,cl) >= 0.6,[1 1 3]));
    h=figure; imshow(id_im);
    %         set(gca,'LooseInset',get(gca,'TightInset'))
    %         set(gcf,'color','white') % White background for the figure.
    %filename = fullfile(mixture_vonMises_dir,[im_splitStr{1},'_cl',num2str(cl),'.png']);
    %imwrite(id_im,filename,'png');
end

x = -pi:0.1:pi;
%c = ['r','g','b'];
c = [ 128 0 128; 205 145 158; 0 0 0]./255;
figure;
%histogram(theta,'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'BinWidth',0.2);
%histogram(theta,'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'NumBins',50);
drawWheel(theta,50,[0.8 0.8 0.8]);
all_mu_hat_polar = [mu_hat_polar, mu_white];
all_kappa_hat = [kappa_hat, kappa_white]; 
hold on;
for cl=1:(numClusters+1)
    yk = prior_probs(cl)*circ_vmpdf(x, all_mu_hat_polar(cl), all_kappa_hat(cl));
    %plot(x, yk,'Color',c(cl),'LineStyle','-','LineWidth',2); hold on;
    circ_line(x,yk,c(cl,:));
end
