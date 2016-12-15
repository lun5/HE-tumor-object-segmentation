% generate output at OIS for best images in the test set
% Luong Nguyen 07/23/2015
% Print outputs based on best region metrics
% UPDATE: 5/31/2016

%% mac
githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir)); cd(githubdir)
seismdir = '/Users/lun5/Research/github/seism'; addpath(genpath(seismdir));
DATA_DIR = '/Users/lun5/Research/data/';
GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','well_defined');
%GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','invasive');

IMG_DIR = fullfile(DATA_DIR,'Tiles_512');%'/home/lun5/HEproject/data/images/test';
%IMG_DIR = fullfile(DATA_DIR,'TilesForLabeling_tiff_renamed');%'/home/lun5/HEproject/data/images/test';
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated','best_images_july30');%fullfile(DATA_DIR,'data','groundTruth_512_512');
%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','eval_reannotated');
%% window
%githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; % window
%seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism';
%addpath(genpath(seismdir)); cd(githubdir)
%DATA_DIR = 'Z:\';
%IMG_DIR = 'Z:\Tiles_512\Test';
%GT_DIR = 'Z:\HEproject\data\groundTruth_512_512';
%RESULTS_DIR = fullfile(DATA_DIR,'HEproject','evaluation_results','Isola_lowres_accurate');
%gt_display = fullfile(DATA_DIR,'groundTruth','groundTruth_reannotated_display');
%outDir = fullfile(RESULTS_DIR,'best_boundary');

%% linux
% githubdir = '/home/lun5/github/HE-tumor-object-segmentation';
% addpath(genpath(githubdir));cd(githubdir);
% seismdir = '/home/lun5/github/seism'; addpath(genpath(seismdir));% linux
% bsrdir = '/home/lun5/github/BSR/grouping';addpath(genpath(bsrdir));
% DATA_DIR ='/home/lun5/HEproject/'; % linux
% IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/Test';
dir_names = {fullfile('eval_non_expert','Maurice'),fullfile('eval_non_expert','Om'),...
    'PMI_lowres_accurate','SuperPixel',fullfile('GraphRLM','new_params'),...
    'GlandSeg','bsr','Isola_speedy',fullfile('JSEG','new_params','scale1'),...
    'ncut_multiscale_1_6','MeanShift',fullfile('EGB','seism_params')};
method_names = {'non-expert','non-expert','colorStats','inNucDist','GraphRLM','GlandSeg','gPb','crisp-bound',...
    'JSEG','NCut','MeanShift','EGB'};

dir_names = {'PMI_lowres_accurate',fullfile('GraphRLM','new_params'),...
    'bsr','Isola_speedy',fullfile('JSEG','new_params','scale1'),...
    'ncut_multiscale_1_6','MeanShift',fullfile('EGB','seism_params')};
method_names = {'colorStats','GraphRLM','gPb','crisp-bound',...
    'JSEG','NCut','MeanShift','EGB'};
RESULTS_DIR = cell(length(dir_names),1);

for i = 1: length(RESULTS_DIR);
    RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',dir_names{i});
end

img_list = dirrec(GT_DIR,'.mat');

ff_scores = zeros(length(RESULTS_DIR),1);
bb_scores = zeros(length(RESULTS_DIR),1);
for med = 1:length(RESULTS_DIR)
    fprintf('\n\nCalculate best segmentation for methods %s...\n',method_names{med}); T = tic;
    %EV_DIR = fullfile(RESULTS_DIR{med},'ev_txt_invasive');
    EV_DIR = fullfile(RESULTS_DIR{med},'ev_txt_invasive_overlap_April4');
    %EV_DIR = fullfile(RESULTS_DIR{med},'ev_txt_well_defined_overlap_April4');
    %eval_bdry_img = dlmread(fullfile(EV_DIR,'eval_bdry_img.txt'));
    %eval_bdry = dlmread(fullfile(EV_DIR,'eval_bdry.txt'));
    load(fullfile(EV_DIR,'fr_mat.mat'));
    fr_mean = mean(fr_mat,1);
    [~,thres_indx] = max(fr_mean);
    UCM_DIR = fullfile(RESULTS_DIR{med},'ucm2');
    SEG_DIR = fullfile(RESULTS_DIR{med},'segmented_images');
    if exist(UCM_DIR,'dir')
        load(fullfile(EV_DIR,'thresh.mat'));
        best_thres = thresh(thres_indx);
    else
        best_thres = thres_indx;
    end
    %% find the optimal partition
    %bdry_outDir = fullfile(RESULTS_DIR{med},'best_bdry_300_May30_well_defined');
    %bdry_outDir = fullfile(RESULTS_DIR{med},'best_bdry_300_May30_invasive');
    %bdry_outDir = fullfile(RESULTS_DIR{med},'best_bdry_May31_overlap_invasive');
    bdry_outDir = fullfile(RESULTS_DIR{med},'best_bdry_May31_overlap_well_defined');
    if ~exist(bdry_outDir,'dir')
        mkdir(bdry_outDir)
    end
    fb_scores = zeros(size(fr_mat,1),1);
    for i = 1:numel(img_list)
        [~,im_name,~] = fileparts(img_list{i}); im_name = lower(im_name);
        %im_name = '4zjh6oyq06xe';
        bdry_outFile = fullfile(bdry_outDir,[im_name, '.tif']);
        I = imread(fullfile(IMG_DIR,[im_name '.tif']));
        if exist(UCM_DIR,'dir')
            if ~exist(fullfile(UCM_DIR,[im_name '.mat']),'file')
                continue;
            end
            tmp = load(fullfile(UCM_DIR,[im_name '.mat']));
            ucm2 = tmp.data;
            labels2 = bwlabel(ucm2 <= best_thres);
            partition = labels2(2:2:end,2:2:end);
            ucm2 = ucm2(3:2:end,3:2:end);
            %bdry_thr = thresh(best_bdry_thres);%(i); best_bdry_thres
            bdry_thr = best_thres;
            bdry_edge_map = (ucm2>=bdry_thr);           
            %Fop_thr = best_Fop_thres(i);
            %Fop_edge_map = (ucm2>=Fop_thr);
        else
            if ~exist(fullfile(SEG_DIR,[im_name '.mat']),'file')
                continue;
            end
            tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
            segs = tmp.data;
            bdry_thr = floor(best_thres); %floor(best_bdry_thres(i));
            bdry_edge_map = logical(seg2bdry(segs{bdry_thr},'imageSize'));
            partition = segs{bdry_thr};
            %Fop_thr = floor(best_Fop_thres(i));
            %Fop_edge_map = edge(segs{Fop_thr});
        end
                
        %% load ground truth
        tmp = load(fullfile(GT_DIR,[im_name '.mat']));
        groundTruth = tmp.groundTruth;
        gt = groundTruth{1}.Segmentation;
        % evaluate the results
        fb = eval_segm(partition, gt, 'fb');
        fb_scores(i) = fb(1);
        %% use ff, bb score        
%         [ff_score, bb_score] = evalRegions(groundTruth,partition);
%         %fr = ff_score*bb_score;
%         if bb_score == -1 % no background in gt
%             fr = ff_score;
%         else
%             alpha = 0.75;
%             fr = alpha*ff_score + (1-alpha)*bb_score;
%         end
%         ff_scores(med) = ff_score; bb_scores(med) = bb_score;
%         %fprintf('Method %s ff_score = %.2f bb_score = %.2f fr = %.2f\n',...
%         %    method_names{med}, ff_score, bb_score,fr);
%         %% print out the images
%         bdry_edge_map = imdilate(bdry_edge_map, strel('disk',2));
%         bdry_edge_map_im = I.*uint8(repmat(~bdry_edge_map,[1 1 3]));
%         pad_im = padarray(bdry_edge_map_im,[60, 60],255,'both');
%         pad_im = insertText(pad_im,[200 0],method_names{med},'FontSize',50,....
%             'BoxColor','white', 'BoxOpacity', 0);
%         pad_im = insertText(pad_im,[50 570],...
%            sprintf('Fb = %.2f, Fr = %.2f',fb(1),fr),'FontSize',50,'BoxColor','white', 'BoxOpacity', 0);
%         imwrite(pad_im,bdry_outFile,'Resolution',300);
    end
    
    %% statistics on the scores
    
    % regions
    fr_scores = fr_mat(:,thres_indx);
    t = toc(T); fprintf('Methods fr_mean=%.2f, fr_std=%.2f, fb_mean=%.2f, fb_std=%.2f, done: %1.2f sec\n',...
        mean(fr_scores),std(fr_scores), mean(fb_scores), std(fb_scores), t);
end
disp('Done');

