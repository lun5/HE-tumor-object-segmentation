% Mixture of von Mises distributions
% This one we calculate the rotation matrix separately for each image
% Luong Nguyen 03/03/2015
%function mixtureVM_allTiles_isolate_white
% filter out white with cluster centers from kmeans
% EM for purple and pink

sourcedir = 'Z:\';
%addpath(genpath(pwd));
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
mixture_vonMises_dir = fullfile(sourcedir,'mixture_von_mises','isolate_white_nonfreeze');

if ~exist(mixture_vonMises_dir,'dir')
    mkdir(mixture_vonMises_dir);
    fileattrib(mixture_vonMises_dir,'+w');
end
% options for mixture model
numClusters = 2; % only purple and pink this time
opts_mixture.noise = 1;

%% get the rotation matrix 
% source image
fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths);% 232
load(fullfile(pwd,'DanTrainingData','cluster_centers.mat')); % cluster center to get white 
white_cluster_indx = 6; % position of white cluster in the array
nColors = 11; % total number of clusters to differentiate from white
mu_white = 2.24; kappa_white = 30; % vM mean and concentration of white
cform = makecform('srgb2lab'); % transformation into Lab space/for isolating white
start_time = zeros(numImages,1); end_time = zeros(numImages,1);
parfor j = 1: numImages
    T(j) = tic;
    imname = imagepaths{j}(1:end-4); 
    %imname = 'dYuoYvkb7Q16Q4p';
    raw_image = imread(fullfile(tiles_dir,[imname '.tif']));
    im_rgb = double(raw_image)./255;
    % identify white pixels
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
    [mm,ii] = min(dist_to_clusters,[],2);
    indx_white = ii == white_cluster_indx;
    
    rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
    rotated_coordinates = rotation_matrix.rotation_matrix*X';
    theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:)); %hue
    theta(:,indx_white) = circ_vmrnd(mu_white, kappa_white, sum(indx_white));

    % 1D von Mises mixture model, freeze, no white version 
    X_cart = [cos(theta); sin(theta)]';
    X_cart_nw = X_cart(~indx_white,:); % nw = not white
    %% Call the function
    [ mu_hat_polar_nw,~, kappa_hat_nw,posterior_probs_nw, prior_probs_nw] =...
        moVM_fixWhite(X_cart_nw,numClusters,opts_mixture);
    posterior_probs = zeros(ncols*nrows,4);
    posterior_probs(~indx_white,[1:2 4]) = posterior_probs_nw;
    posterior_probs(indx_white,3) = 1;
    post_images = reshape(posterior_probs,[nrows ncols 4]);
    mu_hat_polar = [mu_hat_polar_nw, mu_white];
    kappa_hat = [kappa_hat_nw, kappa_white];
    num_pixels = ncols*nrows;
    prior_probs_nw = prior_probs_nw/num_pixels *(num_pixels - sum(indx_white));
    prior_probs = [prior_probs_nw(1:2) sum(indx_white)/num_pixels prior_probs_nw(3)];
    
    save_struct = struct('mu_hat_polar',mu_hat_polar,'kappa_hat',kappa_hat,...
       'posterior_probs',posterior_probs,'prior_probs',prior_probs);
    fname = fullfile(mixture_vonMises_dir,[imname,'_stats.mat']);
    parsave(fname, save_struct);
  
    for cl = 1:(numClusters+1+opts_mixture.noise)
        id_im = raw_image.*uint8(repmat(post_images(:,:,cl) >= 0.5,[1 1 3]));
        filename = fullfile(mixture_vonMises_dir,[imname,'_cl',num2str(cl),'.png']);
        %figure; imshow(id_im);
        imwrite(id_im,filename,'png');
    end
    % histogram of von Mises distribution
    x = -pi:0.1:pi;
    c = [ 128 0 128; 205 145 158; 0 0 0]./255;
    figure;
    drawWheel(theta,50,[0.8 0.8 0.8]);
    hold on;
    for cl=1: (numClusters + 1) % added white
        yk = prior_probs(cl)*circ_vmpdf(x, mu_hat_polar(cl), kappa_hat(cl));
        circ_line(x,yk,c(cl,:));
    end        
    hold off; %xlim([-pi pi]);
    set(gcf,'color','white') % White background for the figure.
    filename = fullfile(mixture_vonMises_dir,[imname,'_hist.png']);
    print(gcf, '-dpng', filename);
    t(j) = toc(T(j));
    fprintf('finish with image %s in %.2f seconds\n', imname, t(j));
end
