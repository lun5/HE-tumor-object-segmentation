%% Script to classify the tissue components
% what the cookcoo is a kernel method?
% Luong Nguyen
% June 16, 2015
% 

%% query all the segments that are carcinoma or vessels
%HOMELMCOMPONENTS = 'Z:\HEproject\data\tissue_components_fromSeg_3channels';
HOMELMCOMPONENTS = 'Z:\HEproject\data\tissue_components_features\tissue_components_groundTruth';
features_dir = fullfile(HOMELMCOMPONENTS,'features');
images_dir = fullfile(HOMELMCOMPONENTS,'images');

component_list = {'carcinoma','fat','vessel','stroma','duct','lymphocyte','white'};
num_components = length(component_list);
% for each component

%comp_name = component_list{1};
%% calculate matrix of features
all_features = cell(num_components,1);
for i = 1: num_components
    comp_name = component_list{i};
    listDir = dir(fullfile(features_dir,['*' comp_name '*']));
    fileNames = {listDir.name}';
    numFiles = length(fileNames);
    features_component = cell(numFiles,1);
    
    for j = 1:numFiles
        fileName = fileNames{j};
        feats = load(fullfile(features_dir,fileName));
        features_component{j} = feats.data(1:end);
    end
    features_component = cat(1,features_component{1:end});
    all_features{i} = features_component;
end

%% calculate labels
all_labels = cell(num_components,1);
all_filenames = cell(num_components,1);
parfor i = 1:num_components
    all_labels{i} = repmat({component_list{i}},length(all_features{i}),1);
    comp_name = component_list{i};
    listDir = dir(fullfile(features_dir,['*' comp_name '*']));
    all_filenames{i} = {listDir.name}';
end

%cmp_ind = 1:5;
cmp_ind = [4,5];
all_features_mat = cat(1,all_features{cmp_ind});
labels = cat(1,all_labels{cmp_ind});
fileNames = cat(1,all_filenames{cmp_ind});
% omit the regions whose area < 5000
ind_area = all_features_mat(:,end) > 5000;
% features_3channels = all_features_mat(ind_area,1:end-1); % 3 channels
% features_hue = all_features_mat(ind_area,1:end - 11); 
% features_brightness = all_features_mat(ind_area,end-11:end-6);
% features_saturation = all_features_mat(ind_area,end-6:end-1);
labels = labels(ind_area); fileNames = fileNames(ind_area); 
features_hue = all_features_mat(ind_area,1:end - 1); 
features = features_hue;
% features = features_brightness;
% features = features_saturation;
% features = cat(2,features_hue,features_brightness);
% features = cat(2,features_hue,features_saturation);
% features = features_3channels;
% %% multiclass classification
% % cv partition
% cv = cvpartition(labels,'k',5);
% rng(1);
% t = templateSVM('Standardize',1);
% [~,numerical_labels] = ismember(labels,component_list);
% Md1 = fitcecoc(features,numerical_labels,'Learners',t);
% CVMd1 = crossval(Md1);
% %pool = parpool;
% options = statset('UseParallel',1);
% oosLoss = kfoldLoss(CVMd1,'Options',options);
% oofLabel = kfoldPredict(CVMd1,'Options',options);
% ConfMat = confusionmat(numerical_labels,oofLabel);
% cm = ConfMat./repmat(sum(ConfMat,2),[1 size(ConfMat,1)]);
% plotConfMat(cm,{component_list{cmp_ind}});

% %% compare a few components


%% use t-SNE
tne_dir = 'C:\Users\luong_nguyen\Documents\GitHub\tSNE_matlab';
addpath(genpath(tne_dir));
no_dims = 2; initial_dims = 7; perpexity = 30;
D = squareform(pdist(all_features_mat,'cosine'));
mappedX = tsne(all_features_mat,[],no_dims,initial_dims, perpexity);
%mappedX = tsne_p(D,[],no_dims); 
%% plot
x = mappedX(:,1); y = mappedX(:,2); %z = mappedX(:,3);
%figure;
h = gscatter(x,y,labels,[],'.');
% for each unique group in 'g', set the ZData property appropriately
gu = unique(labels); 
% for k = 1:numel(gu)
%     zdata = z(strcmp(labels,gu{k}));
%     set(h(k), 'ZData', zdata );
% end
% view(3)

ax = get(gca);
xRange = ax.XLim(2) - ax.XLim(1); yRange = ax.YLim(2) - ax.YLim(1);
hold on
for i = 1:5:length(all_features_mat)
    fileName = fileNames{i}(1:end-4);
    img = imread(fullfile(images_dir,[fileName '.tif']));
    %x = [mappedX(i,1)-size(img,1)/2, mappedX(i,1) + size(img,1)/2];
    %y = [mappedX(i,2)-size(img,2)/2, mappedX(i,2) + size(img,2)/2];
    x = mappedX(i,1)- ax.XLim(1) ; y = mappedX(i,2)-ax.YLim(1);
    if x < 0 || y < 0; fprintf('something is wrong at i = %d',i);end;
    axes('unit','normalized','Position',[x./xRange,y./yRange,0.05,0.05],'Visible','off');
    image(0,0,img);axis off;    
    %place_image(img,mappedX(i,:));
end
hold off;set(gca,'FontSize',20);
set(gca,'color',[0 0 0]);
% 
% % create random dataset
% x = randn(200,1);
% y = randn(200,1);
% z = randn(200,1);
% g = [1*ones(50,1); 2*ones(50,1); 3*ones(50,1); 4*ones(50,1); ];
% % call GSCATTER and capture output argument (handles to lines)
% h = gscatter(x, y, g);
% % for each unique group in 'g', set the ZData property appropriately
% gu = unique(g);
% for k = 1:numel(gu)
%       set(h(k), 'ZData', z( g == gu(k) ));
% end
% view(3)

% D = squareform(pdist(all_features_matrix,'cosine'));
% % set diagonal to 0
% D(sub2ind(size(D),1:size(D,1),1:size(D,2))) = 0;
% [Y,eigvals] = cmdscale(D);
% format short g
% [eigvals eigvals/max(abs(eigvals))]
% 
% Dtriu = D(find(tril(ones(size(D,1)),-1)))';
% maxrelerr = max(abs(Dtriu-pdist(Y(:,1:2))))./max(Dtriu)
% 
% hFigure = figure('Position',[100 100 600 600]);
% colors = cat(1,repmat([1 0 0],length(all_features{1}),1),...
%     repmat([0 0 1], length(all_features{2}),1),...
%     repmat([1 1 0], length(all_features{3}),1));
% figure; scatter(Y(:,1),Y(:,2),2,colors);
% plot(Y(:,1),Y(:,2),'.'); hold on;
% 
% for j = 1:10:numFiles
%    fileName = fileNames{j};
%    img = imread(fullfile(images_dir,[fileName(1:end-4) '.tif']));
%    hAxes = axes('Parent',hFigure,'Unit','pixel',...
%    'Position',[Y(j,1), Y(j,2), size(img,1), size(img,2)]);
%    set(hAxes,'Visible','off');
% end