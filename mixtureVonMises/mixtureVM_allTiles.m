% Mixture of von Mises distributions
% This one we calculate the rotation matrix separately for each image
% Luong Nguyen 04/28/2015
function mixtureVM_allTiles
% script to collect data tiles
% run the code in parallel
%pool = gcp;
addpath(genpath(pwd));
sourcedir = 'Z:\';
% folder storing the tiles of interest
% tiles_dir = fullfile(sourcedir,'TilesForLabeling');
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
mixture_vonMises_dir = fullfile(sourcedir,'mixture_von_mises','same_rot_renamed_images_4');

if ~exist(mixture_vonMises_dir,'dir')
    mkdir(mixture_vonMises_dir);
    fileattrib(mixture_vonMises_dir,'+w');
end
% options for mixture model
numClusters = 3;
opts_mixture.noise = 1;
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

fileNames = dir(fullfile(tiles_dir,'*.tif'));
imagepaths = {fileNames.name}';
numImages = length(imagepaths);% 420
parfor j = 11: 21%numImages
        imname = imagepaths{j}; 
        %imname = '9uixINHtjjiS.tif';
        %imname = 'NemcDj9A7SKH.tif';
        %imname = 'jRh62FQ8hUZWlA.tif';
        im_splitStr = regexp(imname,'\.','split');
        raw_image = double(imread(fullfile(tiles_dir,imname)));
        % change this part so that the parfor can run
        %[f_maps] = getFeatures(double(raw_image),scale,which_features,opts_features);
        r = raw_image(:,:,1); %g = raw_image(:,:,2)./255; b = raw_image(:,:,3)./255;
        X = reshape(raw_image,[size(raw_image,1)*size(raw_image,2),3])';
        rotated_coordinates = rotation_matrix*double(X./255);
        %theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
        theta = atan2(rotated_coordinates(3,:),rotated_coordinates(2,:));
        im_theta = reshape(theta,size(r));
        % Start mixture model
        % X = im_theta(:);
        % X_cart = [rotated_coordinates(2,:); rotated_coordinates(3,:)];
        X_cart = [cos(theta); sin(theta)]';
        %% Call the function
        [ mu_hat_polar,~, kappa_hat,posterior_probs, prior_probs] =...
           moVM(X_cart,numClusters,opts_mixture);
        save_struct = struct('mu_hat_polar',mu_hat_polar,'kappa_hat',kappa_hat,...
            'posterior_probs',posterior_probs,'prior_probs',prior_probs);
        fname = fullfile(mixture_vonMises_dir,[im_splitStr{1},'_stats.mat']);
        parsave(fname, save_struct);
        %[indx_membership, centroids_cart] = spkmeans(X_cart,numClusters);
        % membership
        [~, indx_membership] = max(posterior_probs,[],2); % 4 is the uniform noise
        
        for cl = 1:(numClusters+opts_mixture.noise)
            id_cluster = reshape(indx_membership, size(im_theta));
            id_cluster(id_cluster ~=cl) = 0;
            id_cluster(id_cluster ~=0) = 1;
            id_im = uint8(raw_image).*uint8(repmat(id_cluster,1,1,3));
            %h=figure; imshow(id_im);
            %set(gca,'LooseInset',get(gca,'TightInset'))
            %set(gcf,'color','white') % White background for the figure.
            filename = fullfile(mixture_vonMises_dir,[im_splitStr{1},'_cl',num2str(cl),'.png']);
            imwrite(id_im,filename,'png');
        end
        
        x = -pi:0.1:pi;
        c = ['r','g','b'];
        
        figure;
        %histogram(theta,'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'BinWidth',0.2);
        histogram(theta,'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'NumBins',50);
        %drawWheel(theta,50,[0.8 0.8 0.8]);
        hold on;
        for cl=1:numClusters
            yk = prior_probs(cl)*circ_vmpdf(x, mu_hat_polar(cl), kappa_hat(cl));
            plot(x, yk,'Color',c(cl),'LineStyle','-','LineWidth',2); hold on;
            %circ_line(x,yk,c(cl));
        end
        
%         if opts_mixture.noise
%             yk = prior_probs(numClusters + opts_mixture.noise)./(2*pi);
%             %plot(x, yk,'Color','k','LineStyle','-','LineWidth',2); hold on;
%             circ_line(x,yk,'k');
%         end
        
        hold off; %xlim([-pi pi]);
        set(gcf,'color','white') % White background for the figure.
        filename = fullfile(mixture_vonMises_dir,[im_splitStr{1},'_hist.png']);
        print(gcf, '-dpng', filename);
        display(['finish with image ', imname]);
        close all;
    
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

