% write code to summary the results for each of the methods:
% each image has N number of settings
% each setting --> score board
% save the score for each image 
% calculate the mean at each setting
% choose the setting that have the best values
%% Mac
% githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation';
%addpath(genpath(githubdir)); 
% cd(githubdir)
% seismdir = '/Users/lun5/Research/github/seism'; addpath(genpath(seismdir));
% DATA_DIR = '/Users/lun5/Research/HE_Segmentation/';
% GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','all_files');
% IMG_DIR = fullfile(DATA_DIR,'Tiles_512');%'/home/lun5/HEproject/data/images/test';

%% Windows
githubdir ='C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; 
cd(githubdir)
seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism';  addpath(genpath(seismdir));
DATA_DIR = 'D:\Documents\HE_Segmentation\';
GT_DIR = fullfile(DATA_DIR,'data','groundTruth','coarse_fine_GT_512_512','all_files');
%IMG_DIR = fullfile(DATA_DIR,'Tiles_512');%'/home/lun5/HEproject/data/images/test';

%all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
%    'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg','SuperPixel',...
%    fullfile('eval_non_expert','Om'),fullfile('eval_non_expert','Maurice')};

train_fname = fullfile(githubdir,'otherMethods','train_tiles.txt');
train_table = readtable(train_fname,'ReadVariableNames',false, 'Delimiter',',');
train_table.Properties.VariableNames = {'tile_names','wsi_names'};
im_list = train_table.tile_names;

measures = {%
            'fb'  ,... % Precision-recall for boundaries
            'fop' ,... % Precision-recall for objects and parts
            'fr'  ,... % Precision-recall for regions
            'voi' ,... % Variation of information
            'nvoi',... % Normalized variation of information
            'pri' ,... % Probabilistic Rand index
            'sc'  ,'ssc' ,... % Segmentation covering (two directions)
            'dhd' ,'sdhd',... % Directional Hamming distance (two directions)
            'bgm' ,... % Bipartite graph matching
            'vd'  ,... % Van Dongen
            'bce' ,... % Bidirectional consistency error
            'gce' ,... % Global consistency error
            'lce' ,... % Local consistency error
            };
        
test_fname = fullfile(githubdir,'otherMethods','test_tiles.txt');
test_table = readtable(test_fname,'ReadVariableNames',false, 'Delimiter',',');
test_table.Properties.VariableNames = {'tile_names','wsi_names'};
test_im_list = test_table.tile_names;
%%

all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
     'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg','SuperPixel_distanceBased','SuperPixel_featureBased',...
     'SuperPixel_upsized_new'};
%all_methods = {'GraphRLM','GlandSeg'};
%all_methods = {'SuperPixel_featureBased'};

RESULTS_DIR = cell(length(all_methods),1);

for i = 1:length(all_methods)
    if i <= (length(all_methods) - 2)
        RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
    else
        RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
    end
end

%maxDist = 0.0075;

maxDist_vec = [0.0075, 0.01, 0.015, 0.02];
for mm = 1:length(maxDist_vec)
    maxDist = maxDist_vec(mm);
    RESULTS_DIR = cell(length(all_methods),1);
    for i = 1:length(all_methods)
        RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
        %RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
    end
for i = length(all_methods)
   segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   fprintf('Start with method %s\n',all_methods{i});
   if ~exist(param_scan_dir,'dir')
       mkdir(param_scan_dir)
   end
   tic;
   for j = 1:length(im_list)
      % load the mat file
      im_name = lower(im_list{j});
      if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
              continue;
      end
      tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
      segs = tmp.data{1};
      nSegments = length(segs); % segments 2:2:200
      % load ground truth
      tmp = load(fullfile(GT_DIR,[im_name '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      groundTruth = tmp.groundTruth;
      param_scan_results = cell(nSegments,1);
      parfor k = 1:nSegments
          seg = segs{k,1}; 
          if i == length(all_methods)
            seg = seg(1:4:end,1:4:end);
          end
          if min(seg(:)) ~= 1
              diff = min(seg(:)) - 1;
              seg = seg - diff;
          end
          seg = uint16(seg);
          if mm == 1
          curr_results = zeros(1,length(measures)+1);
          for m = 1:length(measures)
             result = eval_segm(seg, gt, measures{m},maxDist);
             curr_results(m) = result(1);
          end
          [ff_score, bb_score] = evalRegions(groundTruth,seg);
          %fr = ff_score*bb_score;         
          if bb_score == -1 % no background in gt
              curr_results(end) = ff_score;
          else
              alpha = 0.75;
              curr_results(end) = alpha*ff_score + (1-alpha)*bb_score;
          end
          param_scan_results{k} = curr_results;
          else
              result = eval_segm(seg, gt, measures{1},maxDist);
              param_scan_results{k} = result(1);
          end
      end
      param_scan_results = cat(1,param_scan_results{:});
      save(fullfile(param_scan_dir,[im_name '.mat']),'param_scan_results');
   end
   fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end
end


%% Non UCM Scan the region based only
%{
for i = 1:length(all_methods)
   fprintf('Start with method %s\n',all_methods{i});
   segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_region');
   if ~exist(param_scan_dir,'dir')
       mkdir(param_scan_dir)
   end
   tic;
   for j = 1:length(im_list)
      % load the mat file
      if exist(fullfile(param_scan_dir,[im_list{j} '.mat']),'file')
              continue;
      end
      tmp = load(fullfile(segmented_image_dir,[im_list{j} '.mat']));
      if i <= length(all_methods) -2
          segs = tmp.data;
      else
          segs = tmp.data{1}; 
      end
      nSegments = length(segs); % segments 2:2:200
      % load ground truth
      tmp = load(fullfile(GT_DIR,[im_list{j} '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      groundTruth = tmp.groundTruth;
      param_scan_results = cell(nSegments,1);
      parfor k = 1:nSegments
          seg = segs{k,1}; 
          
          if i >= length(all_methods)-1
              seg = seg(1:4:end,1:4:end);
          end
          if min(seg(:)) ~= 1
              diff = min(seg(:)) - 1;
              seg = seg - diff;
          end
          seg = uint16(seg);
          [ff_score, bb_score] = evalRegions(groundTruth,seg);
          %fr = ff_score*bb_score;         
          if bb_score == -1 % no background in gt
              fr_new = ff_score;
          else
              alpha = 0.75;
              fr_new = alpha*ff_score + (1-alpha)*bb_score;
          end
          param_scan_results{k} = fr_new;
      end
      param_scan_results = cat(1,param_scan_results{:});
      save(fullfile(param_scan_dir,[im_list{j} '.mat']),'param_scan_results');
   end
   fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end
%}
%{
for md = 1:length(maxDist_vec)
    maxDist = maxDist_vec(md);
    for i = 1:length(all_methods)
        %fprintf('Start with method %s, maxDist = %.4f\n',all_methods{i},maxDist);
        segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images');
        param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
        % load the mat file
%         tmp = load(fullfile(segmented_image_dir,[im_list{1} '.mat']));
%         if i >= length(all_methods)-1
%             segs = tmp.data{1};
%         else
%             segs = tmp.data;
%         end
        tmp = load(fullfile(param_scan_dir,[im_list{1} '.mat'])); 
        nSegments = size(tmp.param_scan_results,1);
        %nSegments = length(segs); % segments 2:2:200
        fprintf('maxDist: %.4f, method: %s, nSegments = %d\n',...
            maxDist, all_methods{i},nSegments);
    end
end


for md = 1%:length(maxDist_vec)
    maxDist = maxDist_vec(md);

for i = 1:length(all_methods)
   fprintf('Start with method %s, maxDist = %.4f\n',all_methods{i},maxDist);
   segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   if ~exist(param_scan_dir,'dir')
       mkdir(param_scan_dir)
   end
   tic;
   for j = 1:length(im_list)
      % load the mat file
      if exist(fullfile(param_scan_dir,[im_list{j} '.mat']),'file')
              continue;
      end
      tmp = load(fullfile(segmented_image_dir,[im_list{j} '.mat']));
      if i == length(all_methods)
          segs = tmp.data{1};
      else
          segs = tmp.data; 
      end
      nSegments = length(segs); % segments 2:2:200
      % load ground truth
      tmp = load(fullfile(GT_DIR,[im_list{j} '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      param_scan_results = cell(nSegments,1);
      parfor k = 1:nSegments
          seg = segs{k,1}; 
          
          if i == length(all_methods)
              seg = seg(1:4:end,1:4:end);
          end
          if min(seg(:)) ~= 1
              diff = min(seg(:)) - 1;
              seg = seg - diff;
          end
          
          seg = uint16(seg);
          curr_results = zeros(1,length(measures));
          for m = 1:length(measures)
             result = eval_segm(seg, gt, measures{m},maxDist);
             curr_results(m) = result(1);
          end
          param_scan_results{k} = curr_results;
      end
      param_scan_results = cat(1,param_scan_results{:});
      save(fullfile(param_scan_dir,[im_list{j} '.mat']),'param_scan_results');
   end
   fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end
end

%}
%{
for i = 1:length(all_methods)
   %fprintf('Start with method %s\n',all_methods{i});
   tic;
   segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images');
   param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_region');
   tmp = load(fullfile(param_scan_dir,[im_list{1} '.mat']));
   nSegments = size(tmp.param_scan_results,1);
   mean_measures_thres = zeros(nSegments,1);
   
   for j = 1:length(im_list)
       im_name = im_list{j};
       tmp = load(fullfile(param_scan_dir,[im_name '.mat']));
       mean_measures_thres = mean_measures_thres + tmp.param_scan_results;
   end
   
   % find the best thres
   mean_measures_thres = mean_measures_thres./length(im_list);
   [max_val, max_id] = max(mean_measures_thres,[],1); 
   save(fullfile(param_scan_dir, 'average_measures.mat'),...
       'mean_measures_thres','max_val','max_id');
end
%}
%% summary non UCM2 methods

maxDist_vec = [0.0075, 0.01, 0.015, 0.02];
%{
for md = 1:length(maxDist_vec)   
    maxDist = maxDist_vec(md);
    fprintf('Max dist is %.4f\n',maxDist);
    all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
        'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg',...
        'SuperPixel_distanceBased','SuperPixel_featureBased'};
    RESULTS_DIR = cell(length(all_methods),1);
    for i = 1:length(all_methods)
        if i <= (length(all_methods) - 2)
            RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
        else
            RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
        end
    end
    
    for i = 1:length(all_methods)
        %fprintf('Start with method %s\n',all_methods{i});
        tic;
        segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images');
        param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
%         if exist(fullfile(param_scan_dir, 'test_measures.mat'),'file')
%             continue;
%         end
        tmp = load(fullfile(param_scan_dir, 'test_measures.mat'));
        test_measures = tmp.test_measures;
        if size(test_measures,2) == 16
            continue;
        end
%         im_name = im_list{1};
%         tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
%         if i <= length(all_methods)-1
%             segs = tmp.data{1};
%         else
%             segs = tmp.data;
%         end
%         nSegments = length(segs); % segments 2:2:200
%         mean_measures_thres = zeros(nSegments,length(measures));
%         for j = 1:length(im_list)
%             im_name = im_list{j};
%             tmp = load(fullfile(param_scan_dir,[im_name '.mat']));
%             mean_measures_thres = mean_measures_thres + tmp.param_scan_results;
%         end
%         % find the best thres
%         mean_measures_thres = mean_measures_thres./length(im_list);
%         [max_val, max_id] = max(mean_measures_thres,[],1);
%         save(fullfile(param_scan_dir, 'average_measures.mat'),...
%             'mean_measures_thres','max_val','max_id');
        tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));
        max_val = tmp.max_val;
        max_id = tmp.max_id;
        nSegments = size(tmp.mean_measures_thres,1);
%         test_measures = zeros(length(test_im_list),length(measures));
%         for j = 1:length(test_im_list)
%             %load the mat file
%             im_name = test_im_list{j};
%             tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
%             if i == length(all_methods)
%                 segs = tmp.data{1};
%             else
%                 segs = tmp.data;
%             end
%             nSegments = length(segs); % segments 2:2:200
%             %load ground truth
%             tmp = load(fullfile(GT_DIR,[im_name '.mat']));
%             gt = tmp.groundTruth{1}.Segmentation;
%             seg = segs{max_id(1),1};
%             if size(seg,1) > 512
%                 seg = seg(1:4:end,1:4:end);
%             end
%             
%             if min(seg(:)) ~= 1
%                 diff = min(seg(:)) - 1;
%                 seg = seg - diff;
%             end
%             
%             seg = uint16(seg);
%             curr_results = zeros(1,length(measures));
%             
%             for m = 1:length(measures)
%                 result = eval_segm(seg, gt, measures{m},maxDist);
%                 curr_results(m) = result(1);
%             end
%             test_measures(j,:) = curr_results;
%         end
        fr_new = zeros(length(test_im_list),1);
        for j = 1:length(test_im_list)
            im_name = test_im_list{j};
            tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
            if i >= length(all_methods)-1
                segs = tmp.data{1};
            else
                segs = tmp.data;
            end
            nSegments = length(segs); % segments 2:2:200
            %load ground truth
            tmp = load(fullfile(GT_DIR,[im_name '.mat']));
            gt = tmp.groundTruth{1}.Segmentation;
            groundTruth = tmp.groundTruth;
            seg = segs{max_id(1),1};
            if size(seg,1) > 512
                seg = seg(1:4:end,1:4:end);
            end
            
            if min(seg(:)) ~= 1
                diff = min(seg(:)) - 1;
                seg = seg - diff;
            end
            
            seg = uint16(seg);
            [ff_score, bb_score] = evalRegions(groundTruth,seg);
            %fr = ff_score*bb_score;
            if bb_score == -1 % no background in gt
                fr_new(j) = ff_score;
            else
                alpha = 0.75;
                fr_new(j) = alpha*ff_score + (1-alpha)*bb_score;
            end
        end
        test_measures = cat(2,test_measures,fr_new);
        mean_test_measures = mean(test_measures,1);
        save(fullfile(param_scan_dir, 'test_measures.mat'),...
            'test_measures','mean_test_measures');
        %fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
    end
    
    for i = 1:length(all_methods)
        param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
        %fprintf('Start with method %s\n',all_methods{i});
        tmp = load(fullfile(param_scan_dir, 'test_measures.mat'));
        mean_test_measures = tmp.mean_test_measures;
        tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));
        max_val = tmp.max_val;
        max_id = tmp.max_id;
        mean_measure_thres = tmp.mean_measures_thres;
        tmp = load(fullfile(RESULTS_DIR{i},'param_scan_region','average_measures.mat'));
        mean_measures_region = tmp.mean_measures_thres;
        %fprintf('Setting %d of %d, Train results: %.4f, %.4f, %.4f\n',max_id(1), size(tmp.mean_measures_thres,1),max_val(1:3))
        %fprintf('Test results: %.4f, %.4f, %.4f,%.4f\n\n',mean_test_measures(1:3),mean_test_measures(end))
        fprintf('%.4f, %.4f, %.4f,%.4f, %.4f, %.4f,%.4f,%.4f, %d\n',...
            mean_measure_thres(max_id(1),1:3),mean_measures_region(max_id(1)),...
            mean_test_measures(1:3),mean_test_measures(end),max_id(1));
    end
end
%}

%% UCM2 methods region based
%{
%all_methods = {'gPb', 'Isola_speedy',fullfile('colorStats_param_scan','exp_2.5_sig_7')};%,'PMI_lowres_accurate','Isola_lowres_accurate'};
joint_exponent_vec = [1.25, 2,2.5,3];
sig_vec = [0.25, 1, 3, 5, 7];
all_methods = cell(length(joint_exponent_vec)*length(sig_vec),1);
count = 0;

for i = 1:length(joint_exponent_vec)    
    for j = 1:length(sig_vec)
        count = count + 1;
        all_methods{count} = fullfile('colorStats_param_scan', ...
            ['exp' '_' num2str(joint_exponent_vec(i))...
            '_sig_' num2str(sig_vec(j))]);
    end
end

%all_methods = {fullfile('colorStats_param_scan','exp_2.5_sig_7')};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
    RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
end


thres = 0.01:0.02:0.99;

for i = 1:length(all_methods)
   %fprintf('Start with method %s\n',all_methods{i});
   tic;
   param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_region');
   tmp = load(fullfile(param_scan_dir,[im_list{1} '.mat']));
   nSegments = size(tmp.param_scan_results,1);
   mean_measures_thres = zeros(nSegments,1);
   
   for j = 1:length(im_list)
       im_name = im_list{j};
       tmp = load(fullfile(param_scan_dir,[im_name '.mat']));
       mean_measures_thres = mean_measures_thres + tmp.param_scan_results;
   end
   
   % find the best thres
   mean_measures_thres = mean_measures_thres./length(im_list);
   [max_val, max_id] = max(mean_measures_thres,[],1); 
   save(fullfile(param_scan_dir, 'average_measures.mat'),...
       'mean_measures_thres','max_val','max_id');
end

for mm = 1:length(maxDist_vec)
    maxDist = maxDist_vec(mm);
for i = 1:length(all_methods)
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   %fprintf('Start with method %s\n',all_methods{i});
   tmp = load(fullfile(param_scan_dir, 'test_measures.mat'));
   mean_test_measures = tmp.mean_test_measures;
   tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));
   max_val = tmp.max_val; max_id = tmp.max_id;
   mean_measures_thres = tmp.mean_measures_thres;
   tmp = load(fullfile(RESULTS_DIR{i},'param_scan_region','average_measures.mat'));
   mean_measures_region = tmp.mean_measures_thres;
   %fprintf('thres = %.2f, Traing f-scores: %.4f, %.4f, %.4f\n',thres(max_id(1)),max_val(1:3))
   %fprintf('testing f-scores: %.4f, %.4f, %.4f, %.4f\n\n',mean_test_measures(1:3),mean_test_measures(end))
   fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,%.4f, %.2f\n',...
       mean_measures_thres(max_id(1),1:3), mean_measures_region(max_id(1)),...
       mean_test_measures(1:3),mean_test_measures(end),thres(max_id(1)));
end
end
%}

%{
for i = 1:length(all_methods)
    fprintf('Start with method %s.', all_methods{i});
    segmented_image_dir = fullfile(RESULTS_DIR{i},'ucm2');
    param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_region');
    if ~exist(param_scan_dir,'dir')
        mkdir(param_scan_dir)
    end
    tic;
    for j = 1:length(im_list)
        % load the mat file
        im_name = lower(im_list{j});
        if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
            continue;
        end
        tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
        ucm2 = tmp.data;
        ucm2 = ucm2(3:2:end,3:2:end);
        nSegments = length(thres); % segments 2:2:200
        % load ground truth
        tmp = load(fullfile(GT_DIR,[im_name '.mat']));
        gt = tmp.groundTruth{1}.Segmentation;
        groundTruth = tmp.groundTruth;
        param_scan_results = cell(nSegments,1);
        parfor k = 1:nSegments
            seg = ucm2 > thres(k);
            partition = bwlabel(ucm2 <= thres(k));
            %seg = uint16(seg);
            
            [ff_score, bb_score] = evalRegions(groundTruth,partition);
            %fr = ff_score*bb_score;
            if bb_score == -1 % no background in gt
                fr_new = ff_score;
            else
                alpha = 0.75;
                fr_new = alpha*ff_score + (1-alpha)*bb_score;
            end
            param_scan_results{k} = fr_new;
        end
        param_scan_results = cat(1,param_scan_results{:});
        save(fullfile(param_scan_dir,[im_name '.mat']),'param_scan_results');
    end
    fprintf(' Finished in %.2f seconds\n',toc);
end 
%}
%% UCM methods, boundary based
%{
for i = 1:length(all_methods)
   %fprintf('Start with method %s\n',all_methods{i});
   tic;
   segmented_image_dir = fullfile(RESULTS_DIR{i},'ucm2'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_region');
   nSegments = length(thres); % segments 2:2:200
   mean_measures_thres = zeros(nSegments,1);
   
   for j = 1:length(im_list)
       im_name = im_list{j};
       tmp = load(fullfile(param_scan_dir,[im_name '.mat']));
       mean_measures_thres = mean_measures_thres + tmp.param_scan_results;
   end
   
   % find the best thres
   mean_measures_thres = mean_measures_thres./length(im_list);
   [max_val, max_id] = max(mean_measures_thres,[],1); 
   save(fullfile(param_scan_dir, 'average_measures.mat'),...
       'mean_measures_thres','max_val','max_id');
end
%}
%{
for i = length(all_methods)
   fprintf('Start with method %s.', all_methods{i});
   segmented_image_dir = fullfile(RESULTS_DIR{i},'ucm2'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   if ~exist(param_scan_dir,'dir')
       mkdir(param_scan_dir)
   end
   tic;
   for j = 1:length(im_list)
      % load the mat file
      im_name = lower(im_list{j});
      if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
              continue;
      end
      tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
      ucm2 = tmp.data; 
      ucm2 = ucm2(3:2:end,3:2:end);
      nSegments = length(thres); % segments 2:2:200
      % load ground truth
      tmp = load(fullfile(GT_DIR,[im_name '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      param_scan_results = cell(nSegments,1);
      parfor k = 1:nSegments
          seg = ucm2 > thres(k); 
          %seg = uint16(seg);
          curr_results = zeros(1,length(measures));
          for m = 1:length(measures)
             result = eval_segm(seg, gt, measures{m},maxDist);
             curr_results(m) = result(1);
          end
          param_scan_results{k} = curr_results;
      end
      param_scan_results = cat(1,param_scan_results{:});
      save(fullfile(param_scan_dir,[im_name '.mat']),'param_scan_results');
   end
   fprintf(' Finished in %.2f seconds\n',toc);
end
%}

%% UCM summary evaluation results
%all_methods = {'gPb', 'Isola_speedy'};%,'PMI_lowres_accurate'};%,'Isola_lowres_accurate'};
%{
joint_exponent_vec = [1.25, 2,2.5,3];
sig_vec = [0.25, 1, 3, 5, 7];
all_methods = cell(length(joint_exponent_vec)*length(sig_vec),1);
count = 0;

for i = 1:length(joint_exponent_vec)    
    for j = 1:length(sig_vec)
        count = count + 1;
        all_methods{count} = fullfile('colorStats_param_scan', ...
            ['exp' '_' num2str(joint_exponent_vec(i))...
            '_sig_' num2str(sig_vec(j))]);
    end
end

%all_methods = {'gPb', 'Isola_speedy'};%,'PMI_lowres_accurate'};%,'Isola_lowres_accurate'};
all_methods = {fullfile('colorStats_param_scan_2048','exp_2.5_sig_24')};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
% 	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end
thres = 0.01:0.02:0.99;

maxDist_vec = [0.0075, 0.01, 0.015, 0.02];

for mm = 1:length(maxDist_vec)
    maxDist = maxDist_vec(mm);
    RESULTS_DIR = cell(length(all_methods),1);
    if mm > 1
        measures = measures(1);
    end
    for i = 1:length(all_methods)
        RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
        %RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
    end
for i = length(all_methods)
   segmented_image_dir = fullfile(RESULTS_DIR{i},'ucm2'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   fprintf('Start with method %s\n',all_methods{i});
   if ~exist(param_scan_dir,'dir')
       mkdir(param_scan_dir)
   end
   tic;
   for j = 1:length(im_list)
      % load the mat file
      im_name = lower(im_list{j});
      if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
              continue;
      end
      tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
      ucm2 = tmp.data;
      ucm2 = ucm2(1:4:end,1:4:end);
      ucm2 = ucm2(3:2:end,3:2:end);
      nSegments = length(thres); % segments 2:2:200
      % load ground truth
      tmp = load(fullfile(GT_DIR,[im_name '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      groundTruth = tmp.groundTruth;
      param_scan_results = cell(nSegments,1);
      parfor k = 1:nSegments          
          seg = ucm2 > thres(k);
          seg = uint16(seg);
          partition = bwlabel(ucm2 <= thres(k));
          %load ground truth
          tmp = load(fullfile(GT_DIR,[im_name '.mat']));
          gt = tmp.groundTruth{1}.Segmentation;
          groundTruth = tmp.groundTruth;
          if mm == 1
          curr_results = zeros(1,length(measures)+1);
          for m = 1:length(measures)
             result = eval_segm(seg, gt, measures{m},maxDist);
             curr_results(m) = result(1);
          end
          [ff_score, bb_score] = evalRegions(groundTruth,partition);
          %fr = ff_score*bb_score;         
          if bb_score == -1 % no background in gt
              curr_results(end) = ff_score;
          else
              alpha = 0.75;
              curr_results(end) = alpha*ff_score + (1-alpha)*bb_score;
          end
          param_scan_results{k} = curr_results;
          else
              result = eval_segm(seg, gt, measures{1},maxDist);
              param_scan_results{k} = result(1);
          end
      end
      param_scan_results = cat(1,param_scan_results{:});
      save(fullfile(param_scan_dir,[im_name '.mat']),'param_scan_results');
   end
   fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end
end
%}
%{
%maxDist = 0.02;
for mm = 1:length(maxDist_vec)   
    maxDist = maxDist_vec(mm);
    fprintf('Max dist is %.4f\n',maxDist);
    
for i = 1:length(all_methods)
   %fprintf('Start with method %s\n',all_methods{i});
   tic;
   segmented_image_dir = fullfile(RESULTS_DIR{i},'ucm2'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   if exist(fullfile(param_scan_dir, 'test_measures.mat'),'file')
       continue;
   end
   nSegments = length(thres); % segments 2:2:200
   mean_measures_thres = zeros(nSegments,length(measures));
   
   for j = 1:length(im_list)
       im_name = im_list{j};
       tmp = load(fullfile(param_scan_dir,[im_name '.mat']));
       mean_measures_thres = mean_measures_thres + tmp.param_scan_results;
   end
   
   % find the best thres
   mean_measures_thres = mean_measures_thres./length(im_list);
   [max_val, max_id] = max(mean_measures_thres,[],1); 
   save(fullfile(param_scan_dir, 'average_measures.mat'),...
       'mean_measures_thres','max_val','max_id');
%    tmp = load(fullfile(param_scan_dir, 'test_measures.mat'));
%    test_measures = tmp.test_measures;
%    if size(test_measures,2) == 16
%        test_measures(:,16) = [];
%        mean_test_measures = mean(test_measures,1);
%        save(fullfile(param_scan_dir, 'test_measures.mat'),...
%        'test_measures','mean_test_measures');
%        continue;
%    end
   
%    tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));
%    max_val = tmp.max_val;
%    max_id = tmp.max_id;
   
   test_measures = zeros(length(test_im_list),length(measures)+1);
   
   for j = 1:length(test_im_list)
       %load the mat file
       im_name = test_im_list{j};
       if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
           continue;
       end
       tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
       ucm2 = tmp.data;
       ucm2 = ucm2(3:2:end,3:2:end);
       seg = ucm2 > thres(max_id(1)); 
       seg = uint16(seg);
       partition = bwlabel(ucm2 <= thres(max_id(1)));
       %load ground truth
       tmp = load(fullfile(GT_DIR,[im_name '.mat']));
       gt = tmp.groundTruth{1}.Segmentation;
       groundTruth = tmp.groundTruth;
       
       curr_results = zeros(1,length(measures)+1);
       for m = 1:length(measures)
           result = eval_segm(seg, gt, measures{m},maxDist);
           curr_results(m) = result(1);
       end
       % evalRegion
       [ff_score, bb_score] = evalRegions(groundTruth,partition);
       %fr = ff_score*bb_score;
       if bb_score == -1 % no background in gt
           curr_results(end) = ff_score;
       else
           alpha = 0.75;
           curr_results(end)= alpha*ff_score + (1-alpha)*bb_score;
       end       
       test_measures(j,:) = curr_results;
   end
   mean_test_measures = mean(test_measures,1);
   
%    fr_new = zeros(length(test_im_list),1);
%    %test_measures = zeros(length(test_im_list),length(measures));
%    
%    for j = 1:length(test_im_list)
%        %load the mat file
%        %load the mat file
%        im_name = test_im_list{j};
%        if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
%            continue;
%        end
%        tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
%        ucm2 = tmp.data;
%        ucm2 = ucm2(3:2:end,3:2:end);
%        %seg = ucm2 > thres(max_id(1)); 
%        %seg = uint16(seg);
%        seg = bwlabel(ucm2 <= thres(max_id(1)));
%        %load ground truth
%        tmp = load(fullfile(GT_DIR,[im_name '.mat']));
%        gt = tmp.groundTruth{1}.Segmentation;
%        groundTruth = tmp.groundTruth;
%        
%        [ff_score, bb_score] = evalRegions(groundTruth,seg);
%        %fr = ff_score*bb_score;
%        if bb_score == -1 % no background in gt
%            fr_new(j) = ff_score;
%        else
%            alpha = 0.75;
%            fr_new(j) = alpha*ff_score + (1-alpha)*bb_score;
%        end
%    end
%    test_measures = cat(2,test_measures,fr_new);
%    mean_test_measures = mean(test_measures,1);
   
   save(fullfile(param_scan_dir, 'test_measures.mat'),...
       'test_measures','mean_test_measures');
   fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end

for i = 1:length(all_methods)
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   %fprintf('Start with method %s\n',all_methods{i});
   tmp = load(fullfile(param_scan_dir, 'test_measures.mat'));
   mean_test_measures = tmp.mean_test_measures;
   tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));
   max_val = tmp.max_val; max_id = tmp.max_id;
   mean_measures_thres = tmp.mean_measures_thres;
   tmp = load(fullfile(RESULTS_DIR{i},'param_scan_region','average_measures.mat'));
   mean_measures_region = tmp.mean_measures_thres;
   %fprintf('thres = %.2f, Traing f-scores: %.4f, %.4f, %.4f\n',thres(max_id(1)),max_val(1:3))
   %fprintf('testing f-scores: %.4f, %.4f, %.4f, %.4f\n\n',mean_test_measures(1:3),mean_test_measures(end))
   fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,%.4f, %.2f\n',...
       mean_measures_thres(max_id(1),1:3), mean_measures_region(max_id(1)),...
       mean_test_measures(1:3),mean_test_measures(end),thres(max_id(1)));
end

end
%}

%% evaluate HED
%{
all_methods = {'HED'};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	%RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end
maxDist_vec = [0.0075, 0.01, 0.015, 0.02];


for i = length(all_methods)
    fprintf('Start with method %s.', all_methods{i});
    segmented_image_dir = fullfile(RESULTS_DIR{i},'train');
    param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
    if ~exist(param_scan_dir,'dir')
        mkdir(param_scan_dir)
    end
    tic;
    for j = 1:length(im_list)
        % load the mat file
        im_name = lower(im_list{j});
        if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
            continue;
        end
        tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
        E = tmp.fuse; %edges
        nSegments = length(thres); % segments 2:2:200
        % load ground truth
        tmp = load(fullfile(GT_DIR,[im_name '.mat']));
        gt = tmp.groundTruth{1}.Segmentation;
        param_scan_results = cell(nSegments,1);
        parfor k = 1:nSegments
            E1 = double(E>=max(eps,thrs(k)));
            seg = double(bwmorph(E1,'thin',inf));
            seg = imresize(seg(30:end,30:end),[512,512]);
            seg(seg > 0.4) = 1;seg(seg <=0.4) = 0;
            %seg = uint16(seg);
            curr_results = zeros(1,length(measures));
            for m = 1:length(measures)
                result = eval_segm(seg, gt, measures{m},maxDist);
                curr_results(m) = result(1);
            end
            param_scan_results{k} = curr_results;
        end
        param_scan_results = cat(1,param_scan_results{:});
        save(fullfile(param_scan_dir,[im_name '.mat']),'param_scan_results');
    end
    fprintf(' Finished in %.2f seconds\n',toc);
end
%}

thres = 0.01:0.02:0.99;
se = strel('disk', 1);
%{
for mm = 1:length(maxDist_vec)
    maxDist = maxDist_vec(mm);
    for i = length(all_methods)
        fprintf('Start with method %s. MaxDist = .%4f', all_methods{i}, maxDist);
        segmented_image_dir = fullfile(RESULTS_DIR{i},'train');
        param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
        if ~exist(param_scan_dir,'dir')
            mkdir(param_scan_dir)
        end
        tic;
        for j = 1:length(im_list)
            % load the mat file
            im_name = lower(im_list{j});
%             if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
%                 continue;
%             end
            tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
            E = tmp.fuse; %edges
            nSegments = length(thres); % segments 2:2:200
            % load ground truth
            tmp = load(fullfile(GT_DIR,[im_name '.mat']));
            gt = tmp.groundTruth{1}.Segmentation;
            groundTruth = tmp.groundTruth;
            param_scan_results = cell(nSegments,1);
            parfor k = 1:nSegments
                E1 = double(E>=max(eps,thrs(k)));
                E1 = imresize(E1(30:end,30:end),[512,512]);
                bdry = double(bwmorph(E1,'thin',inf));
                bdry = imdilate(bdry, se);
                seg = bwlabel(1-bdry,4);seg(seg == 0) = 1;
                %seg(seg > 0.4) = 1;seg(seg <=0.4) = 0;
                %seg = uint16(seg);
                if mm == 1
                    curr_results = zeros(1,length(measures)+1);
                    for m = 1:length(measures)
                        result = eval_segm(bdry, gt, measures{m},maxDist);
                        curr_results(m) = result(1);
                    end
                    [ff_score, bb_score] = evalRegions(groundTruth,seg);
                    %fr = ff_score*bb_score;
                    if bb_score == -1 % no background in gt
                        curr_results(end) = ff_score;
                    else
                        alpha = 0.75;
                        curr_results(end) = alpha*ff_score + (1-alpha)*bb_score;
                    end    
                else
                    curr_results = eval_segm(bdry, gt, measures{1},maxDist);
                end
                param_scan_results{k} = curr_results;
            end
            param_scan_results = cat(1,param_scan_results{:});
            save(fullfile(param_scan_dir,[im_name '.mat']),'param_scan_results');
        end
        fprintf(' Finished in %.2f seconds\n',toc);
    end
end
%}
%% Put together the results across multiple maxDist

for i = length(all_methods)
    sum_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_summary');
    if ~exist(sum_scan_dir,'dir')
        mkdir(sum_scan_dir);
    end
    for j = 1 :length(im_list)
        im_name = im_list{j};
        if exist(fullfile(sum_scan_dir,[im_name '.mat']),'file')
            continue;
        end
        %nSegments = length(thres);
        nSegments = 40;
        param_scan_results = zeros(nSegments, length(measures) + 4);
        for mm = 1:length(maxDist_vec)
            maxDist = maxDist_vec(mm);
            param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
            tmp = load(fullfile(param_scan_dir,[im_name '.mat']));
            if mm == 1
               param_scan_results(:,1) = tmp.param_scan_results(:,1);
               param_scan_results(:,5:end) = tmp.param_scan_results(:,2:end);
            else
               param_scan_results(:,mm) = tmp.param_scan_results(:,1);
            end
        end
        save(fullfile(sum_scan_dir,[im_name '.mat']),'param_scan_results');
    end
end
%for mm =1:length(maxDist_vec)
%    maxDist = maxDist_vec(mm);

for i = length(all_methods)
   %fprintf('Start with method %s\n',all_methods{i});
   tic;
   %segmented_image_dir = fullfile(RESULTS_DIR{i},'ucm2'); 
   %segmented_image_dir = fullfile(RESULTS_DIR{i},'test');
   segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images');
   %param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_summary');
   %nSegments = length(thres); % segments 2:2:200
   nSegments = 40;
   mean_measures_thres = zeros(nSegments,length(measures)+4);
     
   for j = 1:length(im_list)
       im_name = im_list{j};
       tmp = load(fullfile(param_scan_dir,[im_name '.mat']));
       mean_measures_thres = mean_measures_thres + tmp.param_scan_results;
   end
   % find the best thres
   mean_measures_thres = mean_measures_thres./length(im_list);
   [max_val, max_id] = max(mean_measures_thres,[],1); 
   save(fullfile(param_scan_dir, 'average_measures.mat'),...
       'mean_measures_thres','max_val','max_id');
   for mm = 1:length(maxDist_vec)
       tic;
       maxDist = maxDist_vec(mm);
       %if exist(fullfile(param_scan_dir, ['test_measures_' num2str(maxDist) '.mat']),'file')
       %    continue;
       %end
       
       %if (mm > 1) && (sum(max_id(mm)== max_id(1:mm-1)) > 0)
       %    same_id = find(max_id(mm)== max_id(1:mm-1),1);
       %    copyfile(fullfile(param_scan_dir,['test_measures_' num2str(maxDist_vec(same_id)) '.mat']),...
       %        fullfile(param_scan_dir, ['test_measures_' num2str(maxDist_vec(mm)) '.mat']));
       %    continue;
       %end
       
       test_measures = zeros(length(test_im_list),length(measures)+1);
       for j = 1:length(test_im_list)
           %load the mat file
           %load the mat file
           im_name = test_im_list{j};
           if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
               continue;
           end
           tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
           
           %E = tmp.fuse; %edges
           %E1 = double(E>=max(eps,thres(max_id(mm))));
           %E1 = imresize(E1(30:end,30:end),[512,512]);
           %bdry = double(bwmorph(E1,'thin',inf));
           %bdry = imdilate(bdry, se);
           %seg = bwlabel(1-bdry,4);seg(seg == 0) = 1;
           
           segs = tmp.data{1};
           seg = segs{max_id(mm),1};
           seg = seg(1:4:end,1:4:end);
           seg = uint16(seg);
           %ucm2 = tmp.data;
           %ucm2 = ucm2(1:4:end,1:4:end);
           %ucm2 = ucm2(3:2:end,3:2:end);
           %seg = ucm2 > thres(max_id(1));
           %seg = uint16(seg);
           %partition = bwlabel(ucm2 <= thres(max_id(1)));
           tmp = load(fullfile(GT_DIR,[im_name '.mat']));
           gt = tmp.groundTruth{1}.Segmentation;
           groundTruth = tmp.groundTruth;
           
           curr_results = zeros(1,length(measures)+1);
           for m = 1:length(measures)
               %result = eval_segm(bdry, gt, measures{m},maxDist);
               result = eval_segm(seg, gt, measures{m},maxDist);
               curr_results(m) = result(1);
           end
           [ff_score, bb_score] = evalRegions(groundTruth,seg);
           %[ff_score, bb_score] = evalRegions(groundTruth,partition);
           %fr = ff_score*bb_score;
           if bb_score == -1 % no background in gt
               curr_results(end) = ff_score;
           else
               alpha = 0.75;
               curr_results(end) = alpha*ff_score + (1-alpha)*bb_score;
           end
           test_measures(j,:) = curr_results;
       end
       mean_test_measures = mean(test_measures,1);
       save(fullfile(param_scan_dir, ['test_measures_' num2str(maxDist) '.mat']),...
           'test_measures','mean_test_measures');
       fprintf('Finish with maxDist %.4f in %.2f seconds\n',maxDist,toc);
       %fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
   end
end

%fprintf('MaxDist=%.4f \n',maxDist);
for i = length(all_methods)
   param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_summary');
   fprintf('Start with method %s\n',all_methods{i});
   tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));
   max_val = tmp.max_val; max_id = tmp.max_id;

   tmp = load(fullfile(param_scan_dir,'average_measures.mat'));
   mean_measures_thres = tmp.mean_measures_thres;
   
  for mm =1:length(maxDist_vec)
      maxDist = maxDist_vec(mm);
      curr_max_id = max_id(mm);
      tmp = load(fullfile(param_scan_dir, ['test_measures_' num2str(maxDist) '.mat']));
      mean_test_measures = tmp.mean_test_measures;
      fprintf('%.4f, %.4f, %.4f,%.4f, %.4f, %.4f,%.4f, %.4f,%.2f\n',...
            mean_measures_thres(curr_max_id,[mm,5:6,end]),...
            mean_test_measures(1:3),mean_test_measures(end),thres(curr_max_id));
   end
end
%end
%}

%% summary HED
%{
for mm =1:length(maxDist_vec)
    maxDist = maxDist_vec(mm);
for i = 1:length(all_methods)
   %fprintf('Start with method %s\n',all_methods{i});
   tic;
   %segmented_image_dir = fullfile(RESULTS_DIR{i},'ucm2'); 
   segmented_image_dir = fullfile(RESULTS_DIR{i},'test');
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   if exist(fullfile(param_scan_dir, 'test_measures.mat'),'file')
       continue;
   end
   nSegments = length(thres); % segments 2:2:200
   mean_measures_thres = zeros(nSegments,length(measures));
  
   for j = 1:length(im_list)
       im_name = im_list{j};
       tmp = load(fullfile(param_scan_dir,[im_name '.mat']));
       mean_measures_thres = mean_measures_thres + tmp.param_scan_results;
   end
   % find the best thres
   mean_measures_thres = mean_measures_thres./length(im_list);
   [max_val, max_id] = max(mean_measures_thres,[],1); 
   save(fullfile(param_scan_dir, 'average_measures.mat'),...
       'mean_measures_thres','max_val','max_id');
   test_measures = zeros(length(test_im_list),length(measures));  
   for j = 1:length(test_im_list)
       %load the mat file
       %load the mat file
       im_name = test_im_list{j};
       if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
           continue;
       end
       tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
       %ucm2 = tmp.data;
       %ucm2 = ucm2(3:2:end,3:2:end);
       %seg = ucm2 > thres(max_id(1)); 
       %seg = uint16(seg);
       %load ground truth
       E = tmp.fuse; %edges
       E1 = double(E>=thres(max_id(1)));
       seg = double(bwmorph(E1,'thin',inf));
       seg = imresize(seg(30:end,30:end),[512,512]);
       seg(seg > 0.4) = 1;seg(seg <=0.4) = 0;
       
       tmp = load(fullfile(GT_DIR,[im_name '.mat']));
       gt = tmp.groundTruth{1}.Segmentation;
       curr_results = zeros(1,length(measures));
       for m = 1:length(measures)
           result = eval_segm(seg, gt, measures{m},maxDist);
           curr_results(m) = result(1);
       end
       test_measures(j,:) = curr_results;
   end
   mean_test_measures = mean(test_measures,1);
   save(fullfile(param_scan_dir, 'test_measures.mat'),...
       'test_measures','mean_test_measures');
   %fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);

end

fprintf('MaxDist=%.4f \n',maxDist);
for i = 1:length(all_methods)
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   fprintf('Start with method %s\n',all_methods{i});
   tmp = load(fullfile(param_scan_dir, 'test_measures.mat'));
   mean_test_measures = tmp.mean_test_measures;
   tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));

   max_val = tmp.max_val; max_id = tmp.max_id;
   fprintf('Thres = %.4f, Train results: %.4f, %.4f, %.4f\n',...
       thres(max_id(1)),max_val(1:3))
   fprintf('Test results: %.4f, %.4f, %.4f\n\n',mean_test_measures(1:3))
end
end
%}

