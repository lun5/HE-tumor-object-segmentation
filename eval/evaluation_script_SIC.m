githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir));cd(githubdir);
seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));
DATA_DIR ='/home/lun5/HEproject/';
IMG_DIR = fullfile(DATA_DIR,'data','Tiles_512');

% githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation'; % mac
% seismdir = '/Users/lun5/Research/github/seism'; %mac
% addpath(genpath(seismdir)); cd(githubdir)
% DATA_DIR = '/Users/lun5/Research/data';
% IMG_DIR  = '/Users/lun5/Research/data/Tiles_512';
RESULTS_DIR = fullfile(DATA_DIR,'normalized_evaluation_results','SIC_1_standard');
GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512', 'all_files');
rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
rotation_matrix = rotation_matrix.rotation_matrix;
mu_s = 2; sigma_s = 3;

img_list = dir(fullfile(IMG_DIR,'*.tif'));
img_list = {img_list.name}';
fprintf('Number of images is %d\n',length(img_list));

opts_affinity = setEnvironment_affinity;
opts_affinity.sig = 3;
opts_affinity.features.which_features = {'hue sat'};
which_features = {'hue sat'};
opts_clustering = setEnvironment_clustering;
opts_clustering.display_progress = false;
%mu_s = .1; sigma_s = 1; % for calculating SIC
if (~exist(RESULTS_DIR,'dir'))
    mkdir(RESULTS_DIR);
end
if (~exist(fullfile(RESULTS_DIR,'E_oriented'),'dir'))
    mkdir(fullfile(RESULTS_DIR,'E_oriented'));
end
if (~exist(fullfile(RESULTS_DIR,'ucm2'),'dir'))
    mkdir(fullfile(RESULTS_DIR,'ucm2'));
end
if (~exist(fullfile(RESULTS_DIR,'edgemap'),'dir'))
    mkdir(fullfile(RESULTS_DIR,'edgemap'));
end

parfor i=1:length(img_list)
    [~,im_name,~] = fileparts(img_list{i});im_name = lower(im_name);
    if (~exist(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),'file'))
        fprintf('\n\nCalculate E oriented %s...',im_name); tic;
        I = imread(fullfile(IMG_DIR,[im_name '.tif']));
        im_rgb = double(I);nrows = size(im_rgb,1); ncols = size(im_rgb,2);
        rgb_coords = reshape(im_rgb,[nrows*ncols,3])';       
        [ sic_coords ] = rgb2sic( rgb_coords, mu_s, sigma_s,[]);% rotation_matrix);
        sic1_im = reshape(sic_coords(1,:),[size(im_rgb,1), size(im_rgb,2)]);
        sic2_im = reshape(sic_coords(2,:),[size(im_rgb,1), size(im_rgb,2)]);
        %f_maps = zeros(nrows, ncols, 2);
        f_maps = zeros(nrows, ncols, 1);
	f_maps(:,:,1) = sic1_im;%f_maps(:,:,2) = sic2_im;
        p = learnP_A_B(f_maps,opts_affinity);
        rf = learnPMIPredictor(f_maps,p,[], which_features, opts_affinity);
        Ws = buildW_pmi(f_maps,rf,p,[], which_features, opts_affinity);
        [~,~,~, E_oriented] = graphSegmentation([],{Ws},{size(f_maps(:,:,1))},im_rgb,opts_clustering);
        E_oriented = imresize(E_oriented{1},[nrows ncols]);
        E = max(E_oriented,[],3);
        imwrite(mat2gray(1-E),fullfile(RESULTS_DIR,'edgemap',[im_name,'_edgemap.tif']));
        %imwrite(mat2gray(1-E),fullfile(RESULTS_DIR,'edgemap',[im_name,'_edgemap.tif']),'Resolution',300);
        parsave(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']),E_oriented);
        t = toc; fprintf('done: %1.2f sec\n', t);
    end
end

E_orienteds = cell(1,length(img_list));
parfor i=1:length(img_list)
    [~,im_name,~] = fileparts(img_list{i});im_name = lower(im_name);
    tmp = load(fullfile(RESULTS_DIR,'E_oriented',[im_name '_E_oriented.mat']));
    E_orienteds{i} = tmp.data;
end
%
max_val = max(cellfun(@max_all,E_orienteds));
parfor i=1:length(img_list)
    E_orienteds{i} = E_orienteds{i}/max_val;
end

parfor i=1:length(img_list)
    [~,im_name,~] = fileparts(img_list{i});im_name = lower(im_name);
    if (~exist(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']),'file'))
        fprintf('\n\nCalculate UCM %s...',im_name); T = tic;
        ucm2 = contours2ucm_crisp_boundaries(mat2gray(E_orienteds{i}),'doubleSize');
        parsave(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']),ucm2);
        t = toc(T); fprintf('done: %1.2f sec\n', t);
    end
end

fprintf('Done');
