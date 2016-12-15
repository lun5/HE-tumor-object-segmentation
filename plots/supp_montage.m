%% generate figures for supplementary materials
%% Luong Nguyen 11/9/15
%% window
%{
DATA_DIR = 'Z:\HEproject';
IMG_DIR = 'Z:\Tiles_512\Test';
%GT_DIR = 'Z:\TilesForLabeling_bestImages\bdry_im';
GT_DIR = fullfile(DATA_DIR,'data','GroundTruth','coarse_fine_GT_512_512');
non_expert_dir = 'Z:\HEproject\evaluation_results\eval_non_expert\Om';
SUPP_DIR = 'Z:\HEproject\evaluation_results\supp_materials\eccb_pri_score';
if ~exist(SUPP_DIR,'dir')
    mkdir(SUPP_DIR);
end
%}

%mac
DATA_DIR = '/Users/lun5/Research/data/';
GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512');
IMG_DIR = fullfile(DATA_DIR,'Tiles_512');%'/home/lun5/HEproject/data/images/test';

non_expert_dir = fullfile(DATA_DIR,'evaluation_results','eval_non_expert','Maurice');
%'eval_PMI_hue_offset'
method_dirs = {'PMI_lowres_accurate','SuperPixel',fullfile('GraphRLM','new_params'),...
    'GlandSeg','bsr','Isola_speedy',fullfile('JSEG','new_params','scale1'),...
    'ncut_multiscale_1_6','MeanShift',fullfile('EGB','seism_params')};
method_names = {'color-stats','inter-nuc-dists','GraphRLM','GlandSeg','gPb','crisp-bound',...
    'JSEG','NCut','NCut','MeanShift','EGB'};
RESULTS_DIR = cell(length(method_dirs),1);

for i =1:length(method_dirs)
    RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',method_dirs{i});
end
%pad_dir = 'best_bdry_300_April3_overlap_metric';%'best_bdry_300_April3_fr_fb';   
pad_dir = 'best_bdry_300_May30_well_defined';
%pad_dir = 'best_bdry_300_May30_invasive';
%pad_dir = 'best_bdry_May31_overlap_well_defined';
%pad_dir = 'best_bdry_May31_overlap_invasive';
%pad_dir1 = 'best_bdry_300_May30_invasive';
SUPP_DIR = fullfile(DATA_DIR,'evaluation_results','supp_materials',pad_dir);
if ~exist(SUPP_DIR,'dir')
    mkdir(SUPP_DIR);
end
im_size = [512 512];
img_list = dirrec(fullfile(GT_DIR,'well_defined'),'.mat');
%img_list = dirrec(fullfile(GT_DIR,'invasive'),'.mat');
%img_list = {'13nedzdzfh','bylklqnsvf4d','nfr1icavoojafx','mbdqhorkuxs'};

%img_list = dirrec(fullfile(non_expert_dir,'segmented_images'),'.mat');
parfor i = 1:length(img_list)
    T = tic;
    [~,im_name,~] = fileparts(img_list{i}); 
    images = cell(12,1);im_name = lower(im_name);%(1:end-5));
    outname = fullfile(SUPP_DIR, [im_name '.png']);
%     if exist(outname,'file')
%         continue;
%     end
%     if ~exist(fullfile(GT_DIR,'pad_300_bigger',[im_name '.tif']),'file')
%         continue;
%     end
%     
%     gt_image = imread(fullfile(GT_DIR,'pad_300_bigger',[im_name '.tif']));
    
    if ~exist(fullfile(GT_DIR,'pad_300',[im_name '.tif']),'file')
        continue;
    end
    
    gt_image = imread(fullfile(GT_DIR,'pad_300',[im_name '.tif']));
    if exist(fullfile(non_expert_dir,pad_dir,[im_name '.tif']),'file')
        non_expert_image = imread(fullfile(non_expert_dir,pad_dir,[im_name '.tif']));
    else
        % white image with the non-expert
        non_expert_image = uint8(255*ones(size(gt_image)));
        non_expert_image = insertText(non_expert_image,[200 0],'non-expert','FontSize',50,....
            'BoxColor','white', 'BoxOpacity', 0);
    end
%     he_image = imread(fullfile(IMG_DIR,[im_name '.tif']));
%     pad_im = padarray(imresize(he_image,im_size),[60 60],255,'both'); 
%     pad_im = insertText(pad_im,[150 5],'H&E image','FontSize',50,'BoxColor','white');
%     %imwrite(pad_im,fullfile(IMG_DIR,[im_name '_pad.tif']),'Resolution',300);
%     images{1} = pad_im; %fullfile(IMG_DIR,[im_name '_pad.tif']);
%     pad_im = padarray(imresize(gt_image,im_size),[60 60],255,'both');
%     pad_im = insertText(pad_im,[150 5],'Ground Truth','FontSize',50,'BoxColor','white');
%     %imwrite(pad_im,fullfile(GT_DIR,[im_name '_pad.tif']),'Resolution',300);
%     images{12} = pad_im; %fullfile(GT_DIR,[im_name '_pad.tif']);
    images{2} = non_expert_image;
    images{1} = gt_image;
    for j = 3:12
        %images{j} = imread(fullfile(RESULTS_DIR{j-2},'best_boundary_300_pad_new',[im_name '.tif']));
        if exist(fullfile(RESULTS_DIR{j-2},pad_dir),'dir')
            images{j} = imread(fullfile(RESULTS_DIR{j-2},pad_dir,[im_name '.tif']));
        else
            images{j} = imread(fullfile(RESULTS_DIR{j-2},pad_dir1,[im_name '.tif']));
        end
        %images{j} = imread(fullfile(RESULTS_DIR{j-1},'best_Fop_300_pad',[im_name '.tif']));
    end
    %montage_im = montage(images,'Size',[2 6]);
    %print(fullfile(SUPP_DIR,num2str(i,'%03d_2')),'-dpng','-r300');    
    zImg = cell2mat(reshape(images, [6,2]).');
    %outname = fullfile(SUPP_DIR, [num2str(i,'%03d') '.png']);
    %outname = fullfile(SUPP_DIR, [im_name '.png']);
    imwrite(zImg,outname,'XResolution',150,'YResolution',150);
    %imwrite(zImg,fullfile(SUPP_DIR, [num2str(i) '.tif']),'Resolution',300);
    t = toc(T); fprintf('done with image %s in %1.2f sec\n', im_name,t);
end

% for i = 1:length(img_list)
%     fprintf('\\includegraphics[width=0.95\\textwidth]{supp/best_boundary/%03d.png}\n',i);
% end
% 
% for i = 1:length(img_list)
%     fprintf('\\includegraphics[width=0.95\\textwidth]{supp/best_Fop/%03d.png}\n',i);
% end

% RESULTS_DIR = cell(12,1);
% RESULTS_DIR{1} = fullfile(DATA_DIR,'evaluation_results','eval_PMI_hue_offset');
% RESULTS_DIR{2} = fullfile(DATA_DIR,'evaluation_results','eval_PJoint_hue_fullscale');
% RESULTS_DIR{3} = fullfile(DATA_DIR,'evaluation_results','Isola_lowres_accurate');
% RESULTS_DIR{4} = fullfile(DATA_DIR,'evaluation_results','Isola_speedy');
% RESULTS_DIR{5} = fullfile(DATA_DIR,'evaluation_results','bsr');
% %RESULTS_DIR{6} = fullfile(DATA_DIR,'evaluation_results','JSEG','new_params','scale1');
% RESULTS_DIR{6} = fullfile(DATA_DIR,'evaluation_results','JSEG','new_params','scale2');
% RESULTS_DIR{7} = fullfile(DATA_DIR,'evaluation_results','ncut_multiscale_1_6');
% RESULTS_DIR{8} = fullfile(DATA_DIR,'evaluation_results','EGB','seism_params');
% %RESULTS_DIR{9} = fullfile(DATA_DIR,'evaluation_results','QuadTree');
% RESULTS_DIR{9} = fullfile(DATA_DIR,'evaluation_results','MeanShift');
% RESULTS_DIR{10} = fullfile(DATA_DIR,'evaluation_results','GraphRLM','new_params');
%SUPP_DIR = 'Z:\HEproject\evaluation_results\supp_materials\best_Fop';
% im_size = 512; sp = 50;
% for i = 1:length(img_list)
%     T = tic;
%     [~,im_name,~] = fileparts(img_list{i}); 
%     images = cell(14,1);
%     gt_image = imread(img_list{i});
%     im_name = lower(im_name(1:end-5));
%     he_image = imread(fullfile(IMG_DIR,[im_name '.tif']));
%     images{1} = imresize(he_image,[im_size im_size]);  
%     images{14} = imresize(gt_image,[im_size im_size]);
%     montage_image = zeros((im_size*2+sp), (im_size*7 + sp*6),3);
%     montage_image(1:im_size,1:im_size,:) = images{1};
%     montage_image((im_size+sp+1):(2*im_size+sp),(6*im_size+5*sp+1):(7*im_size+5*sp),:)= images{14};
%     for j = 2:13
%         I = imread(fullfile(RESULTS_DIR{j-1},'best_boundary_300',[im_name '.tif']));
%         images{j} = imresize(I,[im_size im_size]);
%         [c,r] = ind2sub([7 2],j);
%         montage_image((im_size*(r-1) + sp*(r-1)+1):(im_size*r+ sp*(r-1)),(im_size*(c-1)+ sp*(c-1)+1):(im_size*c+ sp*(c-1)),:) = images{j};        
%     end
%     montage(images,'Size',[2 7]);
%     print(h,fullfile(SUPP_DIR,num2str(i,'%03d')),'-dpng','-r300');
%     t = toc(T); fprintf('done with image %s in %1.2f sec\n', im_name,t);
% end

% image_dir = 'Z:\TilesForLabeling_bestImages\bdry_im';
% im_list = dir(fullfile(image_dir,'*.tif'));
% im_list = {im_list.name}';
% images = cell(16,1);
% for i = 8:3:length(im_list)%1:16
%     I = imread(fullfile(image_dir,im_list{i}));
%     images{floor(i/3)-1} = I;
%     %subplot(4,6,i); imshow(I);
% end
% 
% figure; imdisp(images,'Size',[4 6],'Border',0.05);
% print(fullfile('results','annotation_data'),'-dpng','-r300');
% 
% image_dir = 'C:\Users\luong_nguyen\Box Sync\DDIresearch\EBAposter\LuongPoster\images';
% im_list = {'kidney.jpg','pancreas.jpg','lung.jpg'};
% images = cell(3,1);
% for i = 1:3
%     I = imread(fullfile(image_dir,im_list{i}));
%     images{i} = I;
% end
% 
% figure; imdisp(images,'Size',[1 3],'Border',0.02);
% print(fullfile(image_dir,'variation_organs'),'-dpng','-r300');
