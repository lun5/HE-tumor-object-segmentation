%% convert images to superpixels
% javaaddpath(fullfile(pwd,'jardir','myHelloArchive.jar'))
% javaclasspath
% % Call public static functions of the HelloWorld class
% HelloWorld.main('')
% myHelloObject = HelloWorld  % Create an instance of the class through the constructor
% myHelloObject.getVersion()
% myHelloObject.setVersion(-1)
% myHelloObject.getVersion()
% myHelloObject.setVersion(1)
% myHelloObject.getVersion()
% 
% clear java
% javarmpath(fullfile(pwd,'jardir','myHelloArchive.jar'))
% javaaddpath(fullfile(pwd,'jardir','ExtractObjectInformation.jar'));
close all;
addpath(genpath(pwd));
sourcedir = 'Z:\';
mixture_vonMises_dir = fullfile(sourcedir,'mixture_von_mises','same_rot_renamed_images_3');

tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');
superpixel_input_dir = fullfile(sourcedir,'super_pixels_2');

if ~exist(superpixel_input_dir,'dir')
    mkdir(superpixel_input_dir);
    fileattrib(superpixel_input_dir,'+w');
end

numClusters = 3;
fileNames_stats = dir(fullfile(mixture_vonMises_dir,'*.mat'));
imagepaths = {fileNames_stats.name}';
numImages = length(imagepaths);% 420
opts = setEnvironment_affinity;
which_features = opts.features.which_features;
opts_moVM.noise = 1;
imsize = [2048 2048];

parfor j = 1:numImages
        imname = imagepaths{j}; 
        im_splitStr = regexp(imname,'\_','split');
        fname = fullfile(mixture_vonMises_dir,imname);
        stats_mixture = load(fname);
        data = stats_mixture.data;
        posterior_probs = data.posterior_probs;
        [~, indx_membership] = max(posterior_probs,[],2); % 4 is the uniform noise
        indx_membership(indx_membership == 4) = 3; % set noise to white       
        % nuclei purple 1, stroma pink 2, lumen white 3
        label_image = reshape(indx_membership,imsize);
        line_size = [size(label_image) ones(1,size(label_image,2) - 2)*(-1)];
        label_image = cat(1,line_size, label_image);
        fname = fullfile(superpixel_input_dir,[im_splitStr{1},'_k3']);
        dlmwrite(fname,label_image);
        sprintf('Finish with file %s\n',imname) 
end

% ObjectsFileName = fullfile(sourcedir,'super_pixels',...
%     [im_splitStr{1} '_se1_minNuc3_minStr5_minLum10_circle_map']);
% M=[.41 .10 .310 ;  1 0.31 0.8 ; 0.31 0.9 1 ; 0 1 1 ;...
%     1 1 0 ; 1 0 1 ; 1 1 1 ; 0.5 0 0 ; 0 0 0.5 ];
% 
% d1 = dlmread(ObjectsFileName, ',', 1, 0);
% b=label2rgb(d1,M);
% 
% figure, imshow(b);
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'color','white') % White background for the figure.
% 
% figure; imshow(raw_image);
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'color','white') % White background for the figure.
% 
