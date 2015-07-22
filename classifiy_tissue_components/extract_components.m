%% Luong Nguyen
% june 16, 2015
% script to query different types of tissue components

%% query the images
%clearvars; close all;
LABELME = 'C:\Users\luong_nguyen\Documents\GitHub\LabelMeToolbox';
addpath(genpath(LABELME));
data_dir = 'Z:\HEproject\data';
HOMEIMAGES = fullfile(data_dir,'Images'); % you can set here your default folder
HOMEANNOTATIONS = fullfile(data_dir,'Annotations'); % you can set here your default folder
HOMELMSEGMENTS = fullfile(data_dir,'groundTruth_fake');
HOMELMCOMPONENTS = 'Z:\HEproject\data\tissue_components_fromSeg_3channels';
mixture1d_dir = 'Z:\mixture_von_mises\same_rot_renamed_images_FreezeAll';
tiles_dir = 'Z:\TilesForLabeling_tiff_renamed';

%component_list = {'carcinoma','fat','vessel','stroma','duct','white space','lymphocyte'};
%num_components = length(component_list);
if ~exist(HOMELMCOMPONENTS,'dir') 
    mkdir(HOMELMCOMPONENTS);
end

if ~exist(fullfile(HOMELMCOMPONENTS,'images'),'dir') 
    mkdir(fullfile(HOMELMCOMPONENTS,'images'));
end

if ~exist(fullfile(HOMELMCOMPONENTS,'features'),'dir') 
    mkdir(fullfile(HOMELMCOMPONENTS,'features'));
end
%folderlist = { '10feb04_static_cars_techsquare_lot'};
%LMinstall(folderlist, HOMEIMAGES, HOMEANNOTATIONS);

folderlist = {'renamed_images'};
database = LMdatabase(fullfile(HOMEANNOTATIONS,'users','lun5','renamed_images'));

[D,j] = LMquery(database, 'folder','renamed_images');
Nimages = length(D);
opts = setEnvironment_affinity;
which_features = opts.features.which_features;
Nsamples = 10000;numClusters = 3;
opt_cluster = struct('noise',0);
for ndx = 1:Nimages
    T = tic;
    dd = D(ndx);    
    Nobjects = length(dd.annotation.object);
    [~, seg, names] = LM2segments(dd, [], HOMEIMAGES, HOMELMSEGMENTS);
    num_classes = length(names);
    imname = dd.annotation.filename(1:end-4);
    img = imread(fullfile(tiles_dir,[imname '.tif']));
    [nrows, ncols, c] = size(img);
    tmp = load(fullfile(mixture1d_dir,[imname '_stats.mat']));    
    init_params = struct('theta_hat',tmp.data.mu_hat_polar,'kappa_hat',tmp.data.kappa_hat);
    
    for cl = 1:num_classes
        obj_name = names{cl};
        CC = bwconncomp(seg==cl);
        for cp = 1:length(CC.PixelIdxList)
            % recreate image crop from the mask
            BW = zeros(CC.ImageSize);
            BW(CC.PixelIdxList{cp}) = 1;
            [Y,X] = find(BW);
            polygonRegion = img.*repmat(uint8(BW),[1 1 3]);
            % Image crop:
            crop(1) = max(min(X)-2,1);crop(2) = min(max(X)+2,ncols);
            crop(3) = max(min(Y)-2,1);crop(4) = min(max(Y)+2,nrows);
            crop = round(crop);
            imgCrop = polygonRegion(crop(3):crop(4), crop(1):crop(2), :);
            %figure; imshow(imgCrop)
            imgCrop_name = fullfile(HOMELMCOMPONENTS,'images',[imname '_obj' num2str(cp)...
                '_' obj_name '.tif']);
            imwrite(imgCrop,imgCrop_name,'Resolution',300);
            BW_crop = BW(crop(3):crop(4),crop(1):crop(2));
            I = double(imgCrop);
            %pJoint
            f_maps = getFeatures(I./255,1,which_features,opts);
            %% hue feature
            F = f_maps{1}; F(BW_crop == 0) = [];
            area_mask = numel(F);
            if area_mask < 1000 %|| strcmp(dd.annotation.object(j).occluded ,'yes')
                continue;
            end
            [ mu_hat_polar,~, kappa_hat,~, prior_probs] = moVM([cos(F(:)) sin(F(:))],numClusters,opt_cluster,init_params);
            %init_params = struct('theta_hat',mu_hat_polar,'kappa_hat',kappa_hat,'prior_probs',prior_probs);
            init_params.prior_probs = prior_probs;
            [F,ii,jj] = sampleF(f_maps{1},Nsamples,opts,BW_crop);
            [ params,~, prior_probs] = mixture_of_bivariate_VM(F,9,init_params);
            features = [init_params.prior_probs(1:2), init_params.theta_hat, ...
                init_params.kappa_hat,prior_probs(1:5),params.kappa1(1:6)',....
                params.kappa2(1:6)',params.kappa3(1:6)'];
            %% brightness feature
            brightness = f_maps{2}; brightness(BW_crop == 0) = [];
            features = cat(2,features,[mean(brightness),std(brightness),prctile(brightness,25),...
                prctile(brightness,50),prctile(brightness,75)]);
            %% saturation feature
            saturation = f_maps{3}; saturation(BW_crop == 0) = [];
            features = cat(2,features,[mean(saturation),std(saturation),prctile(saturation,25),...
                prctile(saturation,50),prctile(saturation,75),area_mask]);
            
            %plotPMI_theta;
            % save the parameters
            file_features = fullfile(HOMELMCOMPONENTS,'features',[imname '_obj' num2str(cp)...
                '_' obj_name '.mat']);
            parsave(file_features,features);
        end
    end
    t = toc(T);
    fprintf('\n Done with %s in %d secs\n',imname,t);
end
    
%     %%
%     %img = LMimread(dd, 1, HOMEIMAGES); % Load image
%     imname = dd.annotation.filename(1:end-4);
%     img = imread(fullfile(tiles_dir,[imname '.tif']));
%     [nrows, ncols, c] = size(img);
%     %f_maps_full = getFeatures(double(img),1,which_features,opts);
%     % initialize the parameters for the bivariate von Mises distributions
%     %[F,p1,p2] = sampleF(f_maps_full{1},Nsamples,opts);
%     %numClusters = 3;
%     %[ mu_hat_polar,~, kappa_hat,~, ~] = moVM([cos(F(:,1)) sin(F(:,1))],numClusters);
%     tmp = load(fullfile(mixture1d_dir,[imname '_stats.mat']));    
%     init_params = struct('theta_hat',tmp.data.mu_hat_polar,'kappa_hat',tmp.data.kappa_hat);
%     %init_params.theta_hat = tmp.data.mu_hat_polar;
%     %init_params.kappa_hat = tmp.data.kappa_hat;
%     %init_params.prior_probs = tmp.data.prior_probs;
%     for j = 1:Nobjects
%         [X,Y] = getLMpolygon(dd.annotation.object(j).polygon);       
%         BW = poly2mask(double(X),double(Y), nrows, ncols);
%         polygonRegion = img.*repmat(uint8(BW),[1 1 3]);
%         % Image crop:
%         crop(1) = max(min(X)-2,1);crop(2) = min(max(X)+2,ncols);
%         crop(3) = max(min(Y)-2,1);crop(4) = min(max(Y)+2,nrows);
%         crop = round(crop);
%         imgCrop = polygonRegion(crop(3):crop(4), crop(1):crop(2), :);
%         %figure; imshow(imgCrop)
%         imgCrop_name = fullfile(HOMELMCOMPONENTS,'images',[imname '_obj' num2str(j)...
%             '_' dd.annotation.object(j).name '.tif']);
%         imwrite(imgCrop,imgCrop_name,'Resolution',300);
%         BW_crop = BW(crop(3):crop(4),crop(1):crop(2));
%         I = double(imgCrop);
%         %pJoint
%         f_maps = getFeatures(I,1,which_features,opts);
%         F = f_maps{1}; F(BW_crop == 0) = [];
%         area_mask = numel(F);
%         if area_mask < 1000 || strcmp(dd.annotation.object(j).occluded ,'yes')
%             continue;
%         end
%         [ mu_hat_polar,~, kappa_hat,~, prior_probs] = moVM([cos(F(:)) sin(F(:))],numClusters,opt_cluster,init_params);
%         %init_params = struct('theta_hat',mu_hat_polar,'kappa_hat',kappa_hat,'prior_probs',prior_probs);
%         init_params.prior_probs = prior_probs;
%         [F,ii,jj] = sampleF(f_maps{1},Nsamples,opts,BW_crop);
%         [ params,~, prior_probs] = mixture_of_bivariate_VM(F,9,init_params);
%         features = [init_params.prior_probs(1:2), init_params.theta_hat, ...
%             init_params.kappa_hat,prior_probs(1:5),params.kappa1(1:6)',....
%             params.kappa2(1:6)',params.kappa3(1:6)',area_mask];            
%         %plotPMI_theta;        
%         % save the parameters
%         file_features = fullfile(HOMELMCOMPONENTS,'features',[imname '_obj' num2str(j)...
%             '_' dd.annotation.object(j).name '.mat']);
%         parsave(file_features,features);
%     end
%     fprintf('\n Done with %s\n',imname);
% end