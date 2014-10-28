matfiledir = 'C:\Users\luong_nguyen\Dropbox\DDIresearch\crisp-boundaries\DanTrainingData';
if ~ exist(matfiledir,'dir');
    matfiledir = '/Users/lun5/Research/color_deconvolution/results/140926/';
end
svs_name = 'tp10-867-1';
% svs_name = 'tp10-611';
purple_file = load(fullfile(matfiledir,[svs_name 'training_purple.mat']));
training_data_purple =purple_file.training_data_purple;
pink_file = load(fullfile(matfiledir,[svs_name 'training_pink.mat']));
training_data_pink = pink_file.training_data_pink;

%% get the rotation matrix 
% source image
training_data = [training_data_purple(:,:) training_data_pink(:,1:min(6000,size(training_data_pink,2)))];
[U,~,~] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]';
