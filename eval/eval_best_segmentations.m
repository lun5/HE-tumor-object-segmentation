% generate output at OIS for best images in the test set
% Luong Nguyen 07/23/2015
% this is to test why the performanc isn't very good
% UPDATE: 10/12/15

%% mac
githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation';
addpath(genpath(githubdir)); cd(githubdir)
seismdir = '/Users/lun5/Research/github/seism'; addpath(genpath(seismdir));
DATA_DIR = '/Users/lun5/Research/data/';
%GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','well_defined');
GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','invasive');

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

%IMG_DIR = '/home/lun5/HEproject/data/Tiles_512/';
% githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; % window
% seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism';
% addpath(genpath(seismdir)); cd(githubdir)
% DATA_DIR = 'Z:\HEproject';
% IMG_DIR = 'Z:\Tiles_512';%\Test';
% GT_DIR = 'Z:\HEproject\data\groundTruth_512_512';
% RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','EGB');
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512_reannotated_Oct', 'best_images_july30');
%GT_DIR = fullfile(DATA_DIR,'groundTruth','groundTruth_512_512');
%GT_DIR = 'Z:\TilesForLabeling_bestImages\groundTruth_512_512_reannotated_Oct\best_images_july30';
%RESULTS_DIR = cell(12,1);
%GT_DIR = fullfile(DATA_DIR,'data','GroundTruth','coarse_fine_GT_512_512','invasive');
% dir_names = {'eval_PMI_hue_offset','Isola_speedy',fullfile('JSEG','new_params','scale1'),...
%     'Superpixel','GlandSeg','bsr',fullfile('EGB','seism_params'),'MeanShift',...
%     'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),...
%     fullfile('eval_non_expert','Maurice'),fullfile('eval_non_expert','Om'),...
%     'PMI_lowres_accurate'};
% method_names = {'color-stats','crisp-bound','JSEG','inter-nuc-dist','GlandSeg',...
%     'gPb','EGB','MeanShift','NCut','GraphRLM','non-expert','non-expert','color-stats'};
dir_names = {fullfile('eval_non_expert','Maurice'),fullfile('eval_non_expert','Om'),...
    'PMI_lowres_accurate','SuperPixel',fullfile('GraphRLM','new_params'),...
    'GlandSeg','bsr','Isola_speedy',fullfile('JSEG','new_params','scale1'),...
    'ncut_multiscale_1_6','MeanShift',fullfile('EGB','seism_params')};
method_names = {'non-expert','non-expert','colorStats','NNStats','GraphRLM','GlandSeg','gPb','crisp-bound',...
    'JSEG','NCut','MeanShift','EGB'};
RESULTS_DIR = cell(length(dir_names),1);

for i = 1: length(RESULTS_DIR);
    RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',dir_names{i});
end

% RESULTS_DIR{1} = fullfile(DATA_DIR,'evaluation_results','eval_PMI_hue_offset');
% RESULTS_DIR{2} = fullfile(DATA_DIR,'evaluation_results','eval_PJoint_hue_fullscale');
% RESULTS_DIR{3} = fullfile(DATA_DIR,'evaluation_results','Isola_lowres_accurate');
% RESULTS_DIR{4} = fullfile(DATA_DIR,'evaluation_results','Isola_speedy');
% RESULTS_DIR{5} = fullfile(DATA_DIR,'evaluation_results','bsr');
% RESULTS_DIR{6} = fullfile(DATA_DIR,'evaluation_results','JSEG','new_params','scale1');
% RESULTS_DIR{7} = fullfile(DATA_DIR,'evaluation_results','JSEG','new_params','scale2');
% RESULTS_DIR{8} = fullfile(DATA_DIR,'evaluation_results','ncut_multiscale_1_6');
% RESULTS_DIR{9} = fullfile(DATA_DIR,'evaluation_results','EGB','seism_params');
% RESULTS_DIR{10} = fullfile(DATA_DIR,'evaluation_results','QuadTree');
% RESULTS_DIR{11} = fullfile(DATA_DIR,'evaluation_results','MeanShift');
% RESULTS_DIR{12} = fullfile(DATA_DIR,'evaluation_results','GraphRLM','new_params');

%method_names = {'H&E-hue-PMI','H&E-hue-PJoint',...
%    'Isola-lowres-acc','Isola-speedy','gPb',...
%   'JSEG-ss','JSEG-ms','NCut','EGB','QuadTree','MShift','GraphRLM'};

%RESULTS_DIR = fullfile(DATA_DIR,'evaluation_results','Isola_multiscale');
%nSegments = 50;
%thresh = linspace(1/(nSegments+1),1-1/(nSegments+1),nSegments)';
img_list = dirrec(GT_DIR,'.mat');

%img_list = {'13nedzdzfh','bylklqnsvf4d','nfr1icavoojafx','mbdqhorkuxs'};
%set(0,'defaulttextInterpreter','latex');
ff_scores = zeros(length(RESULTS_DIR),1);
bb_scores = zeros(length(RESULTS_DIR),1);
for med = [3 7 8]; %[4 9 12]%1:2%length(RESULTS_DIR)
    fprintf('\n\nCalculate best segmentation for methods %s...\n',method_names{med}); T = tic;
    %EV_DIR = fullfile(RESULTS_DIR{med},'ev_txt_well_defined');
    EV_DIR = fullfile(RESULTS_DIR{med},'ev_txt_invasive');
    %eval_bdry_img = dlmread(fullfile(EV_DIR,'eval_bdry_img.txt'));
    eval_bdry = dlmread(fullfile(EV_DIR,'eval_bdry.txt'));
    %ev_file = fullfile(EV_DIR,'eval_summary_new.txt');
    %if exist(ev_file,'file')
    %    eval_bdry = dlmread(ev_file,'',0,1);
    %else
    %    eval_bdry = dlmread(fullfile(EV_DIR,'eval_summary.txt'),'',0,1);
    %end
    
    %% find the optimal partition
    best_bdry_thres = eval_bdry(1,1);
    best_bdry_F = eval_bdry(:,5);
    %bdry_outDir = fullfile(RESULTS_DIR{med},'best_bdry_300_May30_well_defined');
    bdry_outDir = fullfile(RESULTS_DIR{med},'best_bdry_300_May30_invasive');
    if ~exist(bdry_outDir,'dir')
        mkdir(bdry_outDir)
    end
    %best_bdry_thres = eval_bdry_img(:,2);
    %best_bdry_F = eval_bdry_img(:,4);
    %fprintf('Number of bdry thres is %d\n',length(best_bdry_thres));
    %     eval_Fop_img = dlmread(fullfile(EV_DIR,'eval_Fop_img.txt'));
    %     Fop_outDir = fullfile(RESULTS_DIR{med},'best_Fop_300_pad');
    %     if ~exist(Fop_outDir,'dir')
    %         mkdir(Fop_outDir)
    %     end
    %     best_Fop_thres = eval_Fop_img(:,2);
    %     best_Fop_F = eval_Fop_img(:,4);
    %fprintf('Number of Fop thres is %d\n',length(best_Fop_thres));
    %IMG_EXT = '.tif';
    %img_list = dirrec(IMG_DIR,IMG_EXT);
    UCM_DIR = fullfile(RESULTS_DIR{med},'ucm2');
    SEG_DIR = fullfile(RESULTS_DIR{med},'segmented_images');
    fr_scores = zeros(numel(img_list),1);
    fb_scores = zeros(numel(img_list),1);
    for i = 1:numel(img_list)
        [~,im_name,~] = fileparts(img_list{i}); im_name = lower(im_name);
        %im_name = 'lszomrlgsc5na4q';
        %im_name = 'k1boslrxx7';
        %im_name = '4zjh6oyq06xe';
        bdry_outFile = fullfile(bdry_outDir,[im_name, '.tif']);
        %Fop_outFile = fullfile(Fop_outDir,[im_name, '.tif']);
        %if ~exist(bdry_outFile,'file') || ~exist(Fop_outFile,'file')
        %I = imread(img_list{i});       
        I = imread(fullfile(IMG_DIR,[im_name '.tif']));
        if exist(UCM_DIR,'dir')
            if ~exist(fullfile(UCM_DIR,[im_name '.mat']),'file')
                fr_scores(i) = -1;
                fb_scores(i) = -1;
                continue;
            end
            tmp = load(fullfile(UCM_DIR,[im_name '.mat']));
            ucm2 = tmp.data;
            labels2 = bwlabel(ucm2 <= best_bdry_thres);
            partition = labels2(2:2:end,2:2:end);
            ucm2 = ucm2(3:2:end,3:2:end);
            %bdry_thr = thresh(best_bdry_thres);%(i); best_bdry_thres
            bdry_thr = best_bdry_thres;
            bdry_edge_map = (ucm2>=bdry_thr);           
            %Fop_thr = best_Fop_thres(i);
            %Fop_edge_map = (ucm2>=Fop_thr);
        else
            if ~exist(fullfile(SEG_DIR,[im_name '.mat']),'file')
                fr_scores(i) = -1;
                fb_scores(i) = -1;
                continue;
            end
            tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
            segs = tmp.data;
            bdry_thr = floor(best_bdry_thres); %floor(best_bdry_thres(i));
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
        %[precision,recall, penalty, ri, gce, vi] = evalPrecisionRecall(groundTruth,partition);
        %fr = 2*(precision/penalty * recall/penalty)/...
        %    (precision/penalty + recall/penalty + ((precision+recall)==0));
        %vi = eval_segm(partition, gt, 'voi');
        %% use ff, bb score
        [ff_score, bb_score] = evalRegions(groundTruth,partition);
        %fr = ff_score*bb_score;
        if bb_score == -1 % no background in gt
            fr = ff_score;
        else
            alpha = 0.75;
            fr = alpha*ff_score + (1-alpha)*bb_score;
        end
        ff_scores(med) = ff_score; bb_scores(med) = bb_score;
        fr_scores(i) = fr;
        %fprintf('Method %s ff_score = %.2f bb_score = %.2f fr = %.2f\n',...
        %    method_names{med}, ff_score, bb_score,fr);
        %% print out the images
        bdry_edge_map = imdilate(bdry_edge_map, strel('disk',2));
        bdry_edge_map_im = I.*uint8(repmat(~bdry_edge_map,[1 1 3]));
%         pad_im = padarray(bdry_edge_map_im,[300, 60],255,'both');
%         pad_im = insertText(pad_im,[150 250],method_names{med},'FontSize',50,'BoxColor','white');
%         pad_im = insertText(pad_im,[50 850],...
%             sprintf('P=%.2f,R=%.2f,Pn=%.2f',precision,recall,penalty),'FontSize',50,'BoxColor','white');
%         pad_im = insertText(pad_im,[50 900],...
%            sprintf('Penalized %.2f, %.2f',precision/penalty,recall/penalty),'FontSize',50,'BoxColor','white');
%         pad_im = insertText(pad_im,[50 950],...
%            sprintf('ri = %.2f, gce = %.2f',ri, gce),'FontSize',50,'BoxColor','white');
%         pad_im = insertText(pad_im,[50 1000],...
%            sprintf('Fb = %.2f, vi =%.2f',fb(1),vi),'FontSize',50,'BoxColor','white');
%         pad_im = pad_im(250:end,:,:);
%        pad_im = padarray(bdry_edge_map_im,[60, 60],255,'both');
%        pad_im = insertText(pad_im,[200 0],method_names{med},'FontSize',50,....
%            'BoxColor','white', 'BoxOpacity', 0);
%        pad_im = insertText(pad_im,[50 570],...
%           sprintf('Fb = %.2f, Fr = %.2f',fb(1),fr),'FontSize',50,'BoxColor','white', 'BoxOpacity', 0);
%        imwrite(pad_im,bdry_outFile,'Resolution',300);
%         h = figure; imshow(bdry_edge_map_im);
%         title(method_names{med},'interpreter','none')
%         xlabel(sprintf('$$F_b$$ = %.2f, $$F_r$$ = %.2f',fb(1),fr),'interpreter','latex')
%         set(gca,'FontSize',30);
%         print(h,'test','-dtiff','-r300');
%         close all;

        %             Fop_edge_map = imdilate(Fop_edge_map, strel('disk',1));
        %             Fop_edge_map_im = I.*uint8(repmat(~Fop_edge_map,[1 1 3]));
        %             pad_im = padarray(Fop_edge_map_im,[60 60],255,'both');
        %             pad_im = insertText(pad_im,[150 5],method_names{med},'FontSize',50,'BoxColor','white');
        %             pad_im = insertText(pad_im,[240 570],...
        %                 sprintf('%.2f',best_Fop_F(i)),'FontSize',50,'BoxColor','white');
        %             imwrite(pad_im,Fop_outFile,'Resolution',300);
        %imwrite(label2rgb(labels),fullfile(output_dir,'seg_im',[im_name, '_' num2str(q_thresh), '_seg.jpg']));
        %end
    end
    %t = toc(T); fprintf('done: %1.2f sec\n', t);
    fr_scores(fr_scores == -1) = [];
    fb_scores(fb_scores == -1) = [];
    t = toc(T); fprintf('Methods fr_mean=%.2f, fr_std=%.2f, fb_mean=%.2f, fb_std=%.2f, done: %1.2f sec\n',...
        mean(fr_scores),std(fr_scores), mean(fb_scores), std(fb_scores), t);
end
disp('Done');

% indx_chosen = [1 3:12];%[3 ,4,7, 9,12];
% bb_scores_new = bb_scores(indx_chosen);
% ff_scores_new = ff_scores(indx_chosen);
% methods_new = method_names(indx_chosen);
% colors = distinguishable_colors(length(indx_chosen));
% colors = distinguishable_colors(5);
% shapes = {'o','>','d','s'};
% h = figure; 
% for i = 1:length(bb_scores_new)
% plot(bb_scores_new(i),ff_scores_new(i),'o',...
%    'MarkerEdgeColor','k','MarkerFaceColor',colors(i,:),'MarkerSize',10);
% shape = shapes{floor(i/3) + 1};
% col = colors(mod(i,5)+1,:);
% plot(bb_scores_new(i),ff_scores_new(i),shape,...
%     'MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',20);
% hold on;
% end
% axis([0 1 0 1]);
% legend(methods_new)
% xlabel('$${\mathbf s_b}$$','interpreter','latex','FontSize',50)
% ylabel('$${\mathbf s_f}$$','interpreter','latex','FontSize',50)
% set(gca,'FontSize',25);
% axis square