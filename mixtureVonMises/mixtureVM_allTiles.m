% Mixture of von Mises distributions
% This one we calculate the rotation matrix separately for each image
function mixtureVM_allTiles
% script to collect data tiles
% run the code in parallel
%pool = gcp;
addpath(genpath(pwd));
sourcedir = 'Z:\';
%svs_fnames = dir(fullfile(sourcedir,'svs','*.svs'));
% svs_fnames = dir(fullfile('/Users/lun5/Research/color_deconvolution/aperio_scans','*.svs'));
%svs_fnames = {svs_fnames.name}';
%num_svs = length(svs_fnames);
% folder storing the tiles of interest
% tiles_dir = fullfile(sourcedir,'TilesForLabeling');
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
% training_dir = fullfile(sourcedir,'ColorsTrainingData');
% tiles_dir = '/Users/lun5/Box Sync/TilesForLabeling';
%mixture_vonMises_dir = fullfile(sourcedir,'mixture_von_mises','different_rot_mat');
mixture_vonMises_dir = fullfile(sourcedir,'mixture_von_mises','same_rot_renamed_images_2');

if ~exist(mixture_vonMises_dir,'dir')
    mkdir(mixture_vonMises_dir);
    fileattrib(mixture_vonMises_dir,'+w');
end
%numtiles = 10;
% bubble_svs = {'tp09-16-1.svs','tp09-39-1.svs','tp09-777-1.svs',...
%     'tp09-1813-1.svs','tp10-420-1.svs','tp10-420-2.svs'};

% set environment for feature calculation
% which_features = {'hue opp'};%, 'brightness opp', 'saturation opp'}; 
% opts_features.features.decorrelate = 0;                              % decorrelate feature channels (done separately for each feature type in which_features)?
% opts_features.plot = false;
% scale = 1;

% options for mixture model
numClusters = 3;
%opts_mixture.noise = 1;
matfiledir = fullfile(pwd,'DanTrainingData');
svs_name = 'tp10-867-1';
%svs_name = 'tp10-611';
purple_file = load(fullfile(matfiledir,[svs_name 'training_purple.mat']));
training_data_purple =purple_file.training_data_purple;
pink_file = load(fullfile(matfiledir,[svs_name 'training_pink.mat']));
training_data_pink = pink_file.training_data_pink;

%% get the rotation matrix 
% source image
training_data = [training_data_purple(:,:) training_data_pink(:,1:min(6000,size(training_data_pink,2)))];
[U,~,~] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]';

% parfor i = 1:num_svs
%     svs_fname = svs_fnames{i};
%     if sum(ismember(bubble_svs,svs_fname)) > 0
%         continue;
%     end
%     
%     splitStr = regexp(svs_fname,'\.','split');
%     fileNames = dir(fullfile(tiles_dir,[splitStr{1} '*.' 'tif']));
%     imagepaths = {fileNames.name}';
fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths);% 420
%     %Load training data
%     training_data_purple=load([training_dir filesep splitStr{1} 'training_purple.mat'],'training_data_purple');
%     training_data_purple = training_data_purple.training_data_purple;
%     training_data_pink=load([training_dir filesep splitStr{1} 'training_pink.mat'],'training_data_pink');
%     training_data_pink = training_data_pink.training_data_pink;
%     % calculate the rotation matrix
%     training_data = [training_data_purple(:,1:min(2000,size(training_data_purple,2)))...
%         training_data_pink(:,1:min(8000,size(training_data_pink,2)))];
%     [U,~,~] = svd(training_data,0);
%     rotation_matrix = [-U(:,1) U(:,2:3)]'; 

parfor j = 1:numImages
        imname = imagepaths{j}; 
        im_splitStr = regexp(imname,'\.','split');
        raw_image = double(imread(fullfile(tiles_dir,imname)));
        % change this part so that the parfor can run
        %[f_maps] = getFeatures(double(raw_image),scale,which_features,opts_features);
        r = raw_image(:,:,1)./255; g = raw_image(:,:,2)./255; b = raw_image(:,:,3)./255;
        rotated_coordinates = rotation_matrix*double([r(:)'; g(:)'; b(:)']);
        %theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
        theta = atan2(rotated_coordinates(3,:),rotated_coordinates(2,:));
        im_theta = reshape(theta,size(r));
        % Start mixture model
        % X = im_theta(:);
        % X_cart = [rotated_coordinates(2,:); rotated_coordinates(3,:)];
        X_cart = [cos(theta); sin(theta)]';
        %% Call the function
        numClusters = 3;
        [ mu_hat_polar,mu_hat_cart, kappa_hat,posterior_probs, prior_probs] =...
           moVM(X_cart,numClusters);
        save_struct = struct('mu_hat_polar',mu_hat_polar,'kappa_hat',kappa_hat,...
            'posterior_probs',posterior_probs,'prior_probs',prior_probs);
        fname = fullfile(mixture_vonMises_dir,[im_splitStr{1},'_stats.mat']);
        %parsave(fname, save_struct);
        %[indx_membership, centroids_cart] = spkmeans(X_cart,numClusters);
        % membership
        [~, indx_membership] = max(posterior_probs,[],2); % 4 is the uniform noise
        
        for cl = 1:(numClusters+1)
            id_cluster = reshape(indx_membership, size(im_theta));
            id_cluster(id_cluster ~=cl) = 0;
            id_cluster(id_cluster ~=0) = 1;
            id_im = uint8(raw_image).*uint8(repmat(id_cluster,1,1,3));
            %h=figure; imshow(id_im);
            %set(gcf,'color','white') % White background for the figure.
            filename = fullfile(mixture_vonMises_dir,[im_splitStr{1},'_cl',num2str(cl),'.png']);
            %print(h, '-dpng', filename);
            %imwrite(id_im,filename,'png');
        end
        
        x = -pi:0.1:pi;
        c = ['r','g','b'];
        
        h=figure;
        histogram(theta,'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'BinWidth',0.2);
        hold on;
        for cl=1:numClusters
            yk = prior_probs(cl)*circ_vmpdf(x, mu_hat_polar(cl), kappa_hat(cl));
            plot(x, yk,'Color',c(cl),'LineStyle','-','LineWidth',2); hold on;
        end
        hold off; xlim([-pi pi]);
        set(gcf,'color','white') % White background for the figure.
        filename = fullfile(mixture_vonMises_dir,[im_splitStr{1},'_hist.png']);
        print(h, '-dpng', filename);
        display(['finish with image ', imname]);
        h=[]; close all;
    
end

%%======================================================
%% STEP 1a: Generate data from two 1D distributions.
% clear vars; close all;
% mu1 = 0;      % Mean
% kappa1 = 70;    % kappa
% m1 = 300;      % Number of points
% 
% mu2 = pi;
% kappa2 = 10;
% m2 = 500;
% 
% mu3 = 5*pi/4;
% kappa3 = 0;
% m3 = 100;
% 
% % Generate the data.
% X1 = circ_vmrnd(mu1,kappa1, m1);
% X2 = circ_vmrnd(mu2,kappa2, m2);
% X3 = circ_vmrnd(mu3,kappa3, m3);
% 
% % have to convert this into cartesian coordiates
% X1_cart = [cos(X1) sin(X1)];
% X2_cart = [cos(X2) sin(X2)];
% X3_cart = [cos(X3) sin(X3)];
% 
% X = [X1; X2; X3];
% X_cart = [X1_cart; X2_cart;X3_cart]; 
% d = size(X_cart,2); % dimension
% %%=====================================================
% %% STEP 1b: Plot the data points and their pdfs.
% 
% x = -pi:0.1:pi;
% y1 = circ_vmpdf(x, mu1, kappa1);
% y2 = circ_vmpdf(x, mu2, kappa2);
% y3 = circ_vmpdf(x, mu3, kappa3);
% 
% figure;
% plot(x, y1, 'b-');
% hold on;
% plot(x, y2, 'r-');
% plot(x, y3, 'g-');
% plot(X1, zeros(size(X1)), 'bx', 'markersize', 10);
% plot(X2, zeros(size(X2)), 'rx', 'markersize', 10);
% plot(X3, zeros(size(X3)), 'g+', 'markersize', 10);
% 
% xlim([-pi pi]);
% set(gcf,'color','white') % White background for the figure.
% hold off

