%% script to extract EGB parameters from seism
seism_dir = '/Users/lun5/Research/github/seism'; %mac
addpath(genpath(seism_dir));% linux
egb_par_dir = fullfile(seism_dir, 'results','EGB');
file_list = dir(fullfile(egb_par_dir,'test_bce*'));
file_list = {file_list.name}';
num_par_settings = length(file_list);
% sigma, k, min
params = zeros(num_par_settings,3);

for i = 1:num_par_settings
    [~,fname,~] = fileparts(file_list{i});
    fname_split = strsplit(fname,'_');
    params(i,:) = [str2double(fname_split{3}), str2double(fname_split{4}), str2double(fname_split{5})];
end

save('params_seism.mat','params')