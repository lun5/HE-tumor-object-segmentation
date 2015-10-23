%% script to extract MeanShift parameters from seism
%% Luong Nguyen 10/22/15

seism_dir = '/Users/lun5/Research/github/seism'; %mac
addpath(genpath(seism_dir));% linux
egb_par_dir = fullfile(seism_dir, 'results','MShift');
file_list = dir(fullfile(egb_par_dir,'test_bce*'));
file_list = {file_list.name}';
num_par_settings = length(file_list);
% 'SpatialBandWidth'(sbw), 'RangeBandWidth'(rbw), 'MinimumRegionArea'(mra)
params = zeros(num_par_settings,3);

for i = 1:num_par_settings
    [~,fname,~] = fileparts(file_list{i});
    fname_split = strsplit(fname,'_');
    params(i,:) = [str2double(fname_split{3}), str2double(fname_split{4}), str2double(fname_split{5})];
end

save('params_seism_mshift.mat','params')