% write code to summary the results for each of the methods:
% each image has N number of settings
% each setting --> score board
% save the score for each image 
% calculate the mean at each setting
% choose the setting that have the best values
%% mac
githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation';
%addpath(genpath(githubdir)); 
cd(githubdir)
seismdir = '/Users/lun5/Research/github/seism'; addpath(genpath(seismdir));
DATA_DIR = '/Users/lun5/Research/HE_Segmentation/';
GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','all_files');
IMG_DIR = fullfile(DATA_DIR,'Tiles_512');%'/home/lun5/HEproject/data/images/test';

% all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
%    'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg','SuperPixel_thres',...
%    fullfile('eval_non_expert','Om'),fullfile('eval_non_expert','Maurice')};

%% window

%% linux

%%
maxDist = 0.01; %[0.0075, 0.01, 0.015, 0.02]
maxDist_vec = [0.0075, 0.01, 0.015, 0.02];
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

for i = 1:length(all_methods)
   fprintf('Start with method %s\n',all_methods{i});
   segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   if ~exist(param_scan_dir,'dir')
       mkdir(param_scan_dir)
   end
   tic;
   for j = 1:length(im_list)
      %load the mat file
      im_name = im_list{j};
      if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
              continue;
      end
      tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
      if i == length(all_methods)
          segs = tmp.data{1};
      else
          segs = tmp.data;
      end
      
      nSegments = length(segs); % segments 2:2:200
      %load ground truth
      tmp = load(fullfile(GT_DIR,[im_name '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      param_scan_results = cell(nSegments,1);
      parfor k = 1:nSegments
          seg = segs{k,1}; 
          if size(seg,1) > 512
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
      save(fullfile(param_scan_dir,[im_name '.mat']),'param_scan_results');
   end
   fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end

%% the UCM ones
%{
all_methods = {'gPb', 'Isola_speedy','PMI_lowres_accurate','Isola_lowres_accurate'};

RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end

maxDist = 0.015;

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
 
thres = 0.01:0.02:0.99;
for i = length(all_methods)
   segmented_image_dir = fullfile(RESULTS_DIR{i},'ucm2'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   if ~exist(param_scan_dir,'dir')
       mkdir(param_scan_dir)
   end
   tic;
   for j = 1:length(im_list)
      %load the mat file
      im_name = im_list{j};
      if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
              continue;
      end
      tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
      ucm2 = tmp.data; 
      ucm2 = ucm2(3:2:end,3:2:end);
      nSegments = length(thres); % segments 2:2:200
      %load ground truth
      tmp = load(fullfile(GT_DIR,[im_name '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      param_scan_results = cell(nSegments,1);
      parfor k = 1:nSegments
          seg = ucm2 > thres(k); 
          seg = uint16(seg);
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
   fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end
%}
%% summary of the evaluation data 
%all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
%    'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg','SuperPixel_thres'};

% RESULTS_DIR = cell(length(all_methods),1);
% for i = 1:length(all_methods)
% 	RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
% 	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
% end

maxDist = 0.02; %[0.0075, 0.01, 0.015, 0.02]
test_fname = fullfile(githubdir,'otherMethods','test_tiles.txt');
test_table = readtable(test_fname,'ReadVariableNames',false, 'Delimiter',',');
test_table.Properties.VariableNames = {'tile_names','wsi_names'};
test_im_list = test_table.tile_names;

maxDist_vec = [0.0075, 0.01, 0.015, 0.02];

for md = 1:length(maxDist_vec)
    maxDist = maxDist_vec(md);
%     all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
%         'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg',...
%         'SuperPixel_distanceBased','SuperPixel_featureBased'};
    all_methods = {'gPb', 'Isola_speedy'};
    RESULTS_DIR = cell(length(all_methods),1);
    for i = 1:length(all_methods)
%         if i <= (length(all_methods) - 2)
            RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
%         else
%             RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
%         end
    end

for i = 1:length(all_methods)
   %fprintf('Start with method %s\n',all_methods{i});
   segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
%   if i <= (length(all_methods) - 2)
       fprintf('scp -r %s lun5@idli.csb.pitt.edu:/home/lun5/HEproject/evaluation_results/%s/.\n',...
           fullfile(param_scan_dir, 'average_measures.mat'),...
           fullfile(all_methods{i},['param_scan' '_' num2str(maxDist)]));
       fprintf('scp -r %s lun5@idli.csb.pitt.edu:/home/lun5/HEproject/evaluation_results/%s/.\n',...
           fullfile(param_scan_dir, 'test_measures.mat'),...
           fullfile(all_methods{i},['param_scan' '_' num2str(maxDist)]));
%    else
%        fprintf('scp -r %s lun5@idli.csb.pitt.edu:/home/lun5/HEproject/normalized_evaluation_results/%s/.\n',...
%            fullfile(param_scan_dir, 'average_measures.mat'),...
%            fullfile(all_methods{i},['param_scan' '_' num2str(maxDist)]));
%        fprintf('scp -r %s lun5@idli.csb.pitt.edu:/home/lun5/HEproject/normalized_evaluation_results/%s/.\n',...
%            fullfile(param_scan_dir, 'test_measures.mat'),...
%            fullfile(all_methods{i},['param_scan' '_' num2str(maxDist)]));
%    end
end
end


for md = 1:length(maxDist_vec)
    maxDist = maxDist_vec(md);
    fprintf('Max dist is %.4f\n',maxDist);
    all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
        'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg',...
        'SuperPixel_distanceBased','SuperPixel_featureBased'};

    RESULTS_DIR = cell(length(all_methods),1);
    for i = 1:length(all_methods)
        if i <= (length(all_methods) - 1)
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
   if exist(fullfile(param_scan_dir, 'test_measures.mat'),'file')
       continue;
   end
   im_name = im_list{1};
   tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
   if i == length(all_methods)
       segs = tmp.data{1};
   else
       segs = tmp.data;
   end
   nSegments = length(segs); % segments 2:2:200
   
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
       im_name = test_im_list{j};
       tmp = load(fullfile(segmented_image_dir,[im_name '.mat']));
       if i == length(all_methods)
           segs = tmp.data{1};
       else
           segs = tmp.data;
       end
       
       nSegments = length(segs); % segments 2:2:200
       %load ground truth
       tmp = load(fullfile(GT_DIR,[im_name '.mat']));
       gt = tmp.groundTruth{1}.Segmentation;
       
       seg = segs{max_id(1),1};
       if size(seg,1) > 512
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
       test_measures(j,:) = curr_results;
   end
   mean_test_measures = mean(test_measures,1);
   save(fullfile(param_scan_dir, 'test_measures.mat'),...
       'test_measures','mean_test_measures');
   %fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end

for i = 1:length(all_methods)
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   fprintf('Start with method %s\n',all_methods{i});
   tmp = load(fullfile(param_scan_dir, 'test_measures.mat'));
   mean_test_measures = tmp.mean_test_measures;
   tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));
   max_val = tmp.max_val; 
   fprintf('Train results: %.4f, %.4f, %.4f\n',max_val(1:3))
   fprintf('Test results: %.4f, %.4f, %.4f\n\n',mean_test_measures(1:3))
end

%% summary evaluation results

all_methods = {'gPb', 'Isola_speedy','PMI_lowres_accurate'};%,'Isola_lowres_accurate'};
RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	%RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end

%maxDist = 0.0075;
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
   
   test_measures = zeros(length(test_im_list),length(measures));
   
   for j = 1:length(test_im_list)
       %load the mat file
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
       
       %load ground truth
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