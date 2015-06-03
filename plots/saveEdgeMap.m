%% output edge map
DATA_DIR ='/home/lun5/HEproject/'; % linux
IMG_DIR = '/home/lun5/HEproject/TilesForLabeling_tiff_renamed';
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_hue_512_512_sig3_exp3'); 
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_hue_512_512_sig3_exp2_newAffinity'); 
RESULTS_DIR = '/home/lun5/HEproject/evaluation_results/eval_hue_512_512_sig3_exp2_newAffinity_Merge';
%% read images
IMG_EXT = '.tif';
img_list = dirrec(IMG_DIR,IMG_EXT); thresh = 0.1;
parfor i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        I = imread(img_list{i});
        mult = 4; dz_im = I(1:mult:end,1:mult:end,:);I = double(dz_im);
        %imwrite(I/255,fullfile(RESULTS_DIR,[im_name,'_512_512.tif']),'Resolution',300);
%         if exist(fullfile(RESULTS_DIR,[im_name '_E_oriented.mat']),'file')
%         tmp = load(fullfile(RESULTS_DIR,[im_name '_E_oriented.mat']));
%         E_oriented = tmp.data;
%         E = max(E_oriented,[],3);
%         imwrite(mat2gray(1-E),fullfile(RESULTS_DIR,[im_name '_edgemap.tif']),'Resolution',300);         
%         end
        if exist(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']),'file')
        tmp = load(fullfile(RESULTS_DIR,'ucm2',[im_name '.mat']));
        E_ucm = tmp.data;
        [segmented_image,~] = ucm2colorsegs(E_ucm,I,thresh);
        imwrite(uint8(segmented_image),fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif']),'Resolution',300);         
        end
        fprintf('\n Print segmented image of %s...',im_name); 
end
%  

% parfor i=1:length(img_list)
%         [~,im_name,~] = fileparts(img_list{i});      
%         I = imread(img_list{i});
%         mult = 2; dz_im = I(1:mult:end,1:mult:end,:);I = double(dz_im);
%         imwrite(I/255,fullfile(RESULTS_DIR,'image1024',[im_name,'_512_512.tif']),'Resolution',300);
%         fprintf('\n Print image 1024 of %s...',im_name); 
% end

groundTruthDisp_dir = '/home/lun5/HEproject/groundTruth_imagesc_512';
% montage
for i=1:length(img_list)
        [~,im_name,~] = fileparts(img_list{i});
        filenames = {fullfile(RESULTS_DIR,'image512',[im_name '_512_512.tif']),...
            fullfile(RESULTS_DIR,'segmented_images01',[im_name '_segmentedImage.tif'])};
                 %fullfile(groundTruthDisp_dir,[im_name '_groundtruth.png']),...
                 %fullfile(RESULTS_DIR,'edgemap',[im_name '_edgemap.tif']),...
        if ~exist(filenames{1},'file') || ~exist(filenames{2},'file'); continue; end;
        %filenames = {img_list{i}, fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif'])};
        if ~exist(fullfile(RESULTS_DIR,'montageImageSegment',[im_name '_montage.tif']),'file')
        h = montage(filenames,'Size',[1 2]);
        MyMontage = getframe(gca); %//
        imwrite(MyMontage.cdata,fullfile(RESULTS_DIR,'montageImageSegment01',[im_name '_montage.tif']),'Resolution',300);
        fprintf('\n Print montage image of %s...',im_name); close all;
        end
end