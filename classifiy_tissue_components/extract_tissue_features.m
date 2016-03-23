%% Luong Nguyen
% March 22, 2016
% script to query different types of tissue components

%% query the images
GT_DIR = 'Z:\HEproject\data\GroundTruth\coarse_fine_GT_512_512\well_defined';
imlist = dir(fullfile(GT_DIR,'*.mat'));
imlist = {imlist.name}';
Nimages = length(imlist);
opts = setEnvironment_affinity;
which_features = opts.features.which_features;
Nsamples = 10000;numClusters = 3;
opt_cluster = struct('noise',0);
tiles_dir = 'Z:\HEproject\normalization_512\well_defined';
components_dir = 'Z:\HEproject\data\tissue_components_features\tissue_components_groundTruth';
mixture1d_dir = 'Z:\mixture_von_mises\vM_normalized_images_512';

if ~exist(components_dir,'dir')
    mkdir(components_dir);
end
if ~exist(fullfile(components_dir,'images'),'dir')
    mkdir(fullfile(components_dir,'images'));
end
if ~exist(fullfile(components_dir,'features'),'dir')
    mkdir(fullfile(components_dir,'features'));
end

parfor i = 1:Nimages
    imname = imlist{i}(1:end-4);
    img = imread(fullfile(tiles_dir,[imname '.tif']));
    nrows = size(img,1); ncols = size(img,2);
    fileseg = fullfile(GT_DIR, [imname '.mat']);
    tmp = load(fullfile(GT_DIR,[imname '.mat']));
    gto = tmp.groundTruth{1};
    seg = gto.Segmentation;
    classes = gto.names;
    tmp = load(fullfile(mixture1d_dir,[imname '_stats.mat']));    
    init_params = struct('theta_hat',tmp.data.mu_hat_polar,'kappa_hat',tmp.data.kappa_hat);
    T = tic;
    fprintf('calculate features of image %s ...',imname);
    for cl = 1:length(classes)
        obj_name = gto.names{cl};
        CC = bwconncomp(seg==cl);
        for cp = 1:length(CC.PixelIdxList)
            % recreate image crop from the mask
            BW = zeros(CC.ImageSize);
            BW(CC.PixelIdxList{cp}) = 1;
            area_mask = sum(BW(:));
            if area_mask < 10000 %|| strcmp(dd.annotation.object(j).occluded ,'yes')
                continue;
            end
            [Y,X] = find(BW);
            polygonRegion = img.*repmat(uint8(BW),[1 1 3]);
            % Image crop:
            crop(1) = max(min(X)-2,1);crop(2) = min(max(X)+2,ncols);
            crop(3) = max(min(Y)-2,1);crop(4) = min(max(Y)+2,nrows);
            crop = round(crop);
            imgCrop = polygonRegion(crop(3):crop(4), crop(1):crop(2), :);
            %figure; imshow(imgCrop)
            imgCrop_name = fullfile(components_dir,'images',[imname '_obj' num2str(cp)...
                '_' obj_name '.tif']);
            imwrite(imgCrop,imgCrop_name);
            BW_crop = BW(crop(3):crop(4),crop(1):crop(2));
            I = double(imgCrop);
            %pJoint
            f_maps = getFeatures(I./255,1,which_features,opts);
            %% hue feature
            F = f_maps{1}; F(BW_crop == 0) = [];
            [ mu_hat_polar,~, kappa_hat,~, prior_probs] = moVM([cos(F(:)) sin(F(:))],numClusters,opt_cluster,init_params);
            %init_params = struct('theta_hat',mu_hat_polar,'kappa_hat',kappa_hat,'prior_probs',prior_probs);
            init_params.prior_probs = prior_probs;
            [F,ii,jj] = sampleF(f_maps{1},Nsamples,opts,BW_crop);
            [ params,~, prior_probs] = mixture_of_bivariate_VM(F,9,init_params);
            %features = [init_params.prior_probs(1:2), init_params.theta_hat, ...
            %    init_params.kappa_hat.^(-1),prior_probs(1:5),params.kappa1(1:6)'.^(-1),....
            %    params.kappa2(1:6)'.^(-1),params.kappa3(1:6)'.^(-1)];
            features = [init_params.prior_probs(1:2),prior_probs(1:5)];
            features = [features ,area_mask];
            %% brightness feature
            %brightness = f_maps{2}; brightness(BW_crop == 0) = [];
            %features = cat(2,features,[mean(brightness),std(brightness),prctile(brightness,25),...
            %    prctile(brightness,50),prctile(brightness,75)]);
            %% saturation feature
            %saturation = f_maps{3}; saturation(BW_crop == 0) = [];
            %features = cat(2,features,[mean(saturation),std(saturation),prctile(saturation,25),...
            %    prctile(saturation,50),prctile(saturation,75),area_mask]);
            
            %plotPMI_theta;
            % save the parameters
            file_features = fullfile(components_dir,'features',[imname '_obj' num2str(cp)...
                '_' obj_name '.mat']);
            parsave(file_features,features);
        end
    end
    fprintf(' done in %.2f seconds\n',toc(T));
end
