%% T2 transformation
% Luong Nguyen 8/19/2015
% After multiplication with standard rotation matrix,
% rot 2 and rot 3 are still highly correlated
% using this correlation, can we identify the stain vectors and hence
% deconvolue the new image?

% read the image
github_dir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation';
tiles_dir = fullfile('Z:\','TilesForLabeling_tiff_renamed'); %window
imname = 'dRfMkOErZY.tif';
raw_image = imread(fullfile(tiles_dir, imname));
% sample 10,000 point
num_samples = 1e4;
num_pixels = size(raw_image,1)*size(raw_image,2);
flatten_image = double(reshape(raw_image,[num_pixels 3]))./255;
randindx = randperm(num_pixels,num_samples);
pix_sample = flatten_image(randindx,:)';
% calculate the pca or these points --> specific rotation matrix
%[u,s,v] = svd(pix_sample,'econ');
%custom_rotation_mat = [-u(:,1)'; u(:,2:end)'];
[coeff,score, latent] = pca(flatten_image);
custom_rotation_mat = coeff';
custom_rotated_coordinates = custom_rotation_mat*flatten_image'; %same as score
% examine the plot of rot2 vs rot3 (or sat*cos(h) vs. sat*sin(h))
figure; plot(custom_rotated_coordinates(2,:), custom_rotated_coordinates(3,:),'.');
%%
% multiply the image with standard rotation matrix
tmp = load(fullfile(github_dir,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
rotation_mat = tmp.rotation_matrix;
% plot the training data after transformation
%% if I add some white to the training data?
tmp = load(fullfile(github_dir,'DanTrainingData','tp10-867-1training_pink.mat'));
pink_data = tmp.training_data_pink;
tmp = load(fullfile(github_dir,'DanTrainingData','tp10-867-1training_purple.mat'));
purple_data = tmp.training_data_purple;
training_data = [pink_data(:,1:3:end), purple_data];

training_rotated = rotation_mat*training_data;
sat = sqrt(training_rotated(2,:).^2 + training_rotated(3,:).^2);
% plot rot2 vs. rot3
figure; plot(training_rotated(2,:), training_rotated(3,:),'.');
set(gca,'FontSize',20);xlabel('Rot 2'); ylabel('Rot 3');
%% histogram of the angles
theta = angle(training_rotated(2,:) + 1i*training_rotated(3,:));
figure; histogram(theta,20);
[N,edges,bin] = histcounts(theta);
ind_first_clump = 8;%12;
ind_second_clump = 16;%20;

first_clump = training_data(:,bin <= ind_first_clump);
second_clump = training_data(:, bin > ind_first_clump & bin <= ind_second_clump);
third_clump = training_data(:,bin > ind_second_clump);

figure; scatter3(first_clump(1,:), first_clump(2,:),first_clump(3,:),8,first_clump'./255,'filled');
set(gca,'FontSize',20);xlabel('R'); ylabel('G'); zlabel('B');

figure; scatter3(second_clump(1,:), second_clump(2,:),second_clump(3,:),8,second_clump'./255,'filled');
set(gca,'FontSize',20);xlabel('R'); ylabel('G'); zlabel('B');

figure; scatter3(third_clump(1,:), third_clump(2,:),third_clump(3,:),8,third_clump'./255,'filled');
set(gca,'FontSize',20);xlabel('R'); ylabel('G'); zlabel('B');

for i = 1:9
figure; scatter3(first_clump(1,i:10:end), first_clump(2,i:10:end),first_clump(3,i:10:end),8,first_clump(:,i:10:end)'./255,'filled');
end
for i = 1:9
figure; scatter3(second_clump(1,i:10:end), second_clump(2,i:10:end),second_clump(3,i:10:end),8,second_clump(:,i:10:end)'./255,'filled');
end
%%
first_clump_rotated = training_rotated(:,bin <= ind_first_clump);
second_clump_rotated = training_rotated(:, bin > ind_first_clump & bin <= ind_second_clump);
third_clump_rotated = training_rotated(:, bin > ind_second_clump);

figure; 
plot(second_clump_rotated(2,:), second_clump_rotated(3,:),'.','MarkerEdgeColor',...
    mean(second_clump,2)'./255,'MarkerSize',8,'MarkerFaceColor', mean(second_clump,2)'./255);
axis([-50 100 -60 40]);
hold on;
plot(first_clump_rotated(2,:), first_clump_rotated(3,:),'.','MarkerEdgeColor',...
     mean(first_clump,2)'./255,'MarkerSize',8,'MarkerFaceColor', mean(first_clump,2)'./255);
plot(third_clump_rotated(2,:), third_clump_rotated(3,:),'k.','MarkerSize',8);
set(gca,'FontSize',20);xlabel('Rot 2'); ylabel('Rot 3');
hold off;
legend('clump 1', 'clump 2', 'clump 3');

%%
figure; 
scatter(second_clump_rotated(2,:), second_clump_rotated(3,:),8,second_clump'./255,'filled');
%axis([-50 100 -60 40]);
hold on;
scatter(first_clump_rotated(2,:), first_clump_rotated(3,:),8,first_clump'./255,'filled');
%scatter(third_clump_rotated(2,:), third_clump_rotated(3,:),8,third_clump'./255,'filled');
set(gca,'FontSize',20);xlabel('Rot 2'); ylabel('Rot 3');
hold off;
%legend('clump 1', 'clump 2', 'clump 3');

figure; 
scatter(second_clump_rotated(1,:), second_clump_rotated(2,:),8,second_clump'./255,'filled');
axis([0 450 -40 80 ]);
hold on;
scatter(first_clump_rotated(1,:), first_clump_rotated(2,:),8,first_clump'./255,'filled');
%scatter(third_clump_rotated(2,:), third_clump_rotated(3,:),8,third_clump'./255,'filled');
set(gca,'FontSize',20);xlabel('Rot 1'); ylabel('Rot 2');
hold off;
legend('clump 1', 'clump 2', 'clump 3');



%% if I add some white to the training data?
tmp = load(fullfile(github_dir,'DanTrainingData','tp10-867-1training_pink.mat'));
pink_data = tmp.training_data_pink;
tmp = load(fullfile(github_dir,'DanTrainingData','tp10-867-1training_purple.mat'));
purple_data = tmp.training_data_purple;
training_data = [pink_data(:,1:3:end), purple_data];
%training_data = [pink_data, purple_data];
%training_data_with_white = [training_data,repmat([255;255;255],1,4000)];
%% examine the rotated coordinate: SVD or pca 
% calculate rotation matrix
[coeff,score,latent] = pca(training_data');
rotation_mat = coeff';

%% another way
centered_data = training_data - repmat(mean(training_data,2), [1 size(training_data,2)]);
[u, s, v] = svd(centered_data*centered_data');
rotation_mat = [-v(:,1)'; v(:,2:end)'];
% rotation_mat = coeff';
rotated_coordinates = rotation_mat*flatten_image';
% find out the pca of rot2 vs. rot3
figure; plot(rotated_coordinates(2,:), rotated_coordinates(3,:),'.');
set(gca,'FontSize',20);
% they are actually do not correlate that much 

% from the pca & brightness --> stain vectors
