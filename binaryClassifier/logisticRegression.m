%% script to classify the pairs of pixels
% for each image, 
% count the number of negative example
% take half of that as training data
% take same amount of that for positive examples
% build the logistic model
% what the help is the score
inputDir = 'Z:\HEproject\data\pixelPairs\';
%inputFname = fullfile(inputDir,'luminance_cue.mat');
%inputFname = fullfile(inputDir,'luminance_cue_sig025.mat');
inputFname = fullfile(inputDir,'color_cue.mat');
%inputFname = fullfile(inputDir,'color_cue_sig025.mat');
%inputFname = fullfile(inputDir,'hue_cue.mat');
%inputFname = fullfile(inputDir,'hue_cue_sig025.mat');

load(inputFname,'feat_sample');
numImages = length(feat_sample);
train_set = cell(numImages,1);
test_set = cell(numImages,1);
%collect training data
parfor i = 1:numImages
   data = feat_sample{i};
   negative_ind = find(data(:,1) == 0);
   positive_ind = find(data(:,1) == 1);
   %numSamples = floor(length(negative_ind)/2);
   numSamples = floor(size(data,1)/2);
   perInd = randperm(length(negative_ind),numSamples);
   train_ind = [negative_ind(perInd); positive_ind(perInd)];
   train_set{i} = data(train_ind,:);
   data(train_ind,:) = [];
   test_set{i} = data;
end

train_set = cat(1,train_set{1:end});
test_set = cat(1,test_set{1:end});
% train the model 
fprintf('Hue different model\n Train:');
y = train_set(:,1);
x = -abs(train_set(:,8)); % difference
x = train_set(:,10);
b = glmfit(x,y,'binomial');
yfit = glmval(b,x,'logit');
%plot(-abs(x),y,'o',-abs(x),yfit,'-','LineWidth',2);
%plot(-abs(x),yfit,'o');
yfit_test = glmval(b,abs(test_set(:,4)),'logit');
yfit_label = yfit_test > 0.5;
pr = sum(yfit_label(yfit_label == 1) == test_set(yfit_label == 1,1))/sum(yfit_label);
rec = sum(yfit_label(test_set(:,1) == 1) == test_set(test_set(:,1) == 1,1))/sum(test_set(:,1));
f_score = 2*(pr*rec)/(pr + rec);

[N,edges,bin] = histcounts(x,100);
edges = edges - (edges(2) - edges(1))/2;
xvals = edges(1:end-1)';
ind = unique(bin);
[grpMean,grpStd] = grpstats(yfit,bin,{'mean','std'});
numLevels = ceil(max(log(N)));
cvec = hot(numLevels+4);
figure; scatter(xvals(ind),grpMean,35,cvec(ceil(log(N(ind)))+1,:),'filled');ylim([0 1]);
colormap('hot'); colorbar('Ticks',[0, 1],'TickLabels',{'0',num2str(numLevels+4)});
% hold on;
% plot(xvals(ind), grpMean + 0.02,'r-');
% plot(xvals(ind), grpMean - 0.02,'b-');
% hold off;

yfit_labels = yfit > 0.5;
pr = sum(yfit_labels(yfit_labels==1) == y(yfit_labels==1))/sum(yfit_labels);
rec = sum(yfit_labels(y==1) == y(y==1))/sum(y);
f_score = 2*(pr*rec)/(pr + rec);
% x = [2100 2300 2500 2700 2900 3100 ...
%      3300 3500 3700 3900 4100 4300]';
% y = [1 1 1 1 1 0 1 0 0 0 0 0]';
% b = glmfit(x,y,'binomial');
% yfit = glmval(b,x,'logit');
% plot(x,y,'o',x,yfit,'-','LineWidth',2);
