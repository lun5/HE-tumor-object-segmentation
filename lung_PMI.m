source_dir = 'Z:\Lung';
tiles_dir = fullfile(source_dir,'tiles');
imnames = {'13fsu-18_im71_stx20351_sty8005.jpg','13fsu-18_im84_stx16281_sty10006.jpg'};

numImages = length(imnames);
opts = setEnvironment_affinity;
which_features = opts.features.which_features{1};

for i = 1:numImages
    %tic;
    imname = fullfile(tiles_dir,imnames{i});
    raw_image = imread(imname);
    ndown = 1;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
    %I = double(dz_im);%  figure; imshow(I/255);size(I)
    figure; imshow(raw_image);
    %% construct raw PMI
    Nsamples = 50000; % to construct the raw histogram
    f_maps_full = getFeatures(double(dz_im),1,{which_features},opts);
    % initialize the parameters for the bivariate von Mises distributions
    [F,~,~] = sampleF(f_maps_full{1},Nsamples,opts);
    numClusters = 3;
    [ mu_hat_polar,~, kappa_hat,~, prior_probs] = moVM([cos(F(:,1)) sin(F(:,1))],numClusters);
    init_params.theta_hat = mu_hat_polar;
    init_params.kappa_hat = kappa_hat;
    init_params.prior_probs = prior_probs;
    %% PMI for components
    [nrows, ncols, c] = size(dz_im);
    h = impoly;
    pos = getPosition(h);
    X = pos(:,1); Y = pos(:,2);
    BW = poly2mask(double(X),double(Y), nrows, ncols);
    polygonRegion = dz_im.*repmat(uint8(BW),[1 1 3]);
    
    crop(1) = max(min(X)-2,1);crop(2) = min(max(X)+2,ncols);
    crop(3) = max(min(Y)-2,1);crop(4) = min(max(Y)+2,nrows);
    crop = round(crop);
                
    % Image crop:
    imgCrop = polygonRegion(crop(3):crop(4), crop(1):crop(2), :);
    figure; imshow(imgCrop)
    BW_crop = BW(crop(3):crop(4),crop(1):crop(2));
    I = double(imgCrop);
    
    f_maps = getFeatures(I,1,{which_features},opts);
    Nsamples = 10000;
    [F,~,~] = sampleF(f_maps{1},Nsamples,opts,BW_crop);
    
    % figure; imshow(uint8(BW_crop)*255); hold on;
    % plot(p1(:,2), p2(:,1),'bx','MarkerSize',10);
    % plot(p2(:,2), p2(:,1),'rx','MarkerSize',10);
    % hold off;
    
    [ params,~, prior_probs] = mixture_of_bivariate_VM(F,9,init_params);
    mixture_params.params = params;
    mixture_params.prior_probs = prior_probs;
    mixture_params.init_params = init_params;
    plotPMI_theta;
end
