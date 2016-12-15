% write code to summary the results for each of the methods:
% each image has N number of settings
% each setting --> score board
% save the score for each image 
% calculate the mean at each setting
% choose the setting that have the best values
githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation';
%addpath(genpath(githubdir)); 
cd(githubdir)
seismdir = '/Users/lun5/Research/github/seism'; addpath(genpath(seismdir));
DATA_DIR = '/Users/lun5/Research/HE_Segmentation/';
GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','all_files');
IMG_DIR = fullfile(DATA_DIR,'Tiles_512');%'/home/lun5/HEproject/data/images/test';

%all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
%    'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg','SuperPixel',...
%    fullfile('eval_non_expert','Om'),fullfile('eval_non_expert','Maurice')};

all_methods = {'MeanShift', fullfile('EGB','seism_params'),fullfile('JSEG','new_params','scale1'), ...
    'ncut_multiscale_1_6',fullfile('GraphRLM','new_params'),'GlandSeg'};

RESULTS_DIR = cell(length(all_methods),1);
for i = 1:length(all_methods)
	%RESULTS_DIR{i} = fullfile(DATA_DIR,'normalized_evaluation_results',all_methods{i});
	RESULTS_DIR{i} = fullfile(DATA_DIR,'evaluation_results',all_methods{i});
end

maxDist = 0.0075;

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
   segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images'); 
   param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan');
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
      segs = tmp.data; 
      nSegments = length(segs); % segments 2:2:200
      % load ground truth
      tmp = load(fullfile(GT_DIR,[im_list{j} '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      param_scan_results = cell(nSegments,1);
      parfor k = 1:nSegments
          seg = segs{k,1}; 
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


