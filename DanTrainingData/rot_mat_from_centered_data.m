% this script test a method for generating rotation matrix,
% in which the data is centered to 0 by subtracting the mean
% Luong Nguyen
% 8/22/2015


%% Directories
github_dir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation';
tiles_dir = fullfile('Z:\','TilesForLabeling_tiff_renamed'); %window

%% training data
%% if I add some white to the training data?
tmp = load(fullfile(github_dir,'DanTrainingData','tp10-867-1training_pink.mat'));
pink_data = tmp.training_data_pink;
tmp = load(fullfile(github_dir,'DanTrainingData','tp10-867-1training_purple.mat'));
purple_data = tmp.training_data_purple;
data = [pink_data(:,1:3:end), purple_data]';
% center data
num_dim = size(data,2);
% dc_vec = ones(1,3)*(1/sqrt(num_dim));
% dc_comp = data*dc_vec';
% mean_im = mean(data - dc_comp*dc_vec,1);
mean_im = mean(data,1);
mean_im = mean_im./norm(mean_im);
mean_comp = data*mean_im';
%left_over = 1/sqrt(num_dim)*(data - dc_comp*dc_vec - mean_comp*mean_im);
left_over = data - mean_comp*mean_im;
[u, s, v] = svd(left_over,'econ');
rotation_mat = v';
rotation_mat = cat(1,rotation_mat(1:2,:),abs(rotation_mat(3,:)));
% Note that the last row of v' is the basis for null space of centered_data'(mx3)
% in this case, it is the same as mean_training
% the brightness is hence
tr_brightness = rotation_mat(3,:)*data';
data_rotated = rotation_mat*data';
tr_sat = sqrt(data_rotated(1,:).^2 + data_rotated(2,:).^2);
% plot rot1 vs. rot2
%figure; plot(training_rotated(1,:), training_rotated(2,:),'.');
figure; scatter(data_rotated(1,:), data_rotated(2,:),8,data./255,'filled');
set(gca,'FontSize',20);xlabel('Rot 1'); ylabel('Rot 2');
tr_hue = angle(data_rotated(1,:) + 1i*data_rotated(2,:));
figure; rose(tr_hue);
figure; histogram(tr_hue,20,'Normalization','probability');
xlim([-pi pi]); set(gca,'FontSize',20);xlabel('Hue');
 % shift the distribution by pi/2 to manually separate the clumps
figure; histogram(mod(tr_hue-pi/2,2*pi),21,'Normalization','probability');
xlim([0 2*pi]); set(gca,'FontSize',20);xlabel('Hue'); % 

[N,edges,bin] = histcounts(tr_hue,20);
[N,edges,bin] = histcounts(mod(tr_hue+pi/2,2*pi),20);
ind_first_clump = 10;%12;
ind_second_clump = max(bin);

data = data';
first_clump = data(:,bin <= ind_first_clump);
second_clump = data(:, bin > ind_first_clump & bin <= ind_second_clump);

% first_clump = data(:,hue_nw <= 0);
% second_clump = data(:,hue_nw> 0);
stp = 100;
figure; scatter3(first_clump(1,1:stp:end), first_clump(2,1:stp:end),first_clump(3,1:stp:end),8,first_clump(:,1:stp:end)'./255,'filled');
set(gca,'FontSize',20);xlabel('R'); ylabel('G'); zlabel('B');

figure; scatter3(second_clump(1,1:stp:end), second_clump(2,1:stp:end),second_clump(3,1:stp:end),8,second_clump(:,1:stp:end)'./255,'filled');
set(gca,'FontSize',20);xlabel('R'); ylabel('G'); zlabel('B');

%%
first_clump_rotated = data_rotated(:,bin <= ind_first_clump);
second_clump_rotated = data_rotated(:, bin > ind_first_clump & bin <= ind_second_clump);

figure; 
scatter(first_clump_rotated(1,1:stp:end), first_clump_rotated(2,1:stp:end),8,first_clump(:,1:stp:end)'./255,'filled');
%axis([-100 100 -100 100]);
set(gca,'FontSize',20);xlabel('Rot 1'); ylabel('Rot 2');
hold on;
stp = 100;
scatter(second_clump_rotated(1,1:stp:end), second_clump_rotated(2,1:stp:end),8,second_clump(:,1:stp:end)'./255,'filled');
%axis([-100 100 -100 100]);
set(gca,'FontSize',20);xlabel('Rot 1'); ylabel('Rot 2');
hold off;
%legend('clump 1', 'clump 2', 'clump 3');

figure; 
scatter(data_rotated(1,1:stp:end)./max(tr_sat), data_rotated(2,1:stp:end)./max(tr_sat),8,...
    data(:,1:stp:end)'./255,'filled');
hold on; plot([0 0 ],[-1 1],'b-'); plot([-1 1],[0 0],'b-'); 
plot(cos(-pi:0.01:pi),sin(-pi:0.01:pi),'r.');
hold off; axis square
set(gca,'FontSize',20);xlabel('Rot 1/ max sat'); ylabel('Rot 2/max sat');

%% Image
imname = 'dRfMkOErZY.tif';
imname = 'Os6RmI2IU30i.tif';
raw_image = imread(fullfile(tiles_dir, imname));
figure; imshow(raw_image);
num_pixels = size(raw_image,1)*size(raw_image,2);
flatten_image = double(reshape(raw_image,[num_pixels 3]))./255;
rotated_coordinates = rotation_mat*flatten_image';
brightness = rotated_coordinates(3,:);
sat = sqrt(rotated_coordinates(1,:).^2 + rotated_coordinates(2,:).^2);
hue = angle(rotated_coordinates(1,:) + 1i*rotated_coordinates(2,:));
% index of white pixels
ind_white = rotated_coordinates(3,:) > sum(rotation_mat(3,:)) - eps;
%ind_white = brightness
flatten_nw = flatten_image(~ind_white,:);
data_rotated = rotated_coordinates(:,~ind_white);
brightness_nw = brightness(~ind_white);
sat_nw = sat(~ind_white);
hue_nw = hue(~ind_white);

figure; rose(hue_nw);
figure; histogram(mod(hue_nw-pi/2,2*pi), 50,'Normalization','probability');
xlim([0 2*pi]); xlabel('hue'); set(gca,'FontSize',20);

[N,edges,bin] = histcounts(mod(hue_nw-pi/2,2*pi), 20);
ind_first_clump = 7;%12;
ind_second_clump = max(bin);

%data = flatten_nw'.*255;
data = flatten_image'.*255; data_rotated = rotated_coordinates;
figure; 
stp = 100;
scatter(data_rotated(1,1:stp:end)./max(sat_nw), data_rotated(2,1:stp:end)./max(sat_nw),8,...
    data(:,1:stp:end)'./255,'filled');
ax = axis;
hold on; plot([0 0],[-1 1],[-1 1],[0 0]); axis tight;
plot(cos(-pi:0.01:pi),sin(-pi:0.01:pi),'r.');
hold off; axis square
set(gca,'FontSize',20);%xlabel('Rot 1/ max sat'); ylabel('Rot 2/max sat');

%% another set of training data
tmp = load(fullfile(pwd,'DanTrainingData','tp10-611training_pink.mat'));
pink_data_2 = tmp.training_data_pink;
tmp = load(fullfile(pwd,'DanTrainingData','tp10-611training_purple.mat'));
purple_data_2 = tmp.training_data_purple;
data_2 = [pink_data_2, purple_data_2]';

tr_brightness = rotation_mat(3,:)*data_2'./255;
data_rotated_2 = rotation_mat*data_2'./255;
tr_sat = sqrt(data_rotated_2(1,:).^2 + data_rotated_2(2,:).^2);
% plot rot1 vs. rot2
%figure; plot(training_rotated(1,:), training_rotated(2,:),'.');
figure; scatter(data_rotated_2(1,:), data_rotated_2(2,:),8,data_2./255,'filled');
set(gca,'FontSize',20);xlabel('Rot 1'); ylabel('Rot 2');
tr_hue = angle(data_rotated_2(1,:) + 1i*data_rotated_2(2,:));
figure; rose(tr_hue);
figure; histogram(tr_hue,20,'Normalization','probability');
xlim([-pi pi]); set(gca,'FontSize',20);xlabel('Hue');

% what if I try to center in the rot 1 & rot 2 plane