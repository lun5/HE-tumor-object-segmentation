%% script to extract parameters for evaluation from seism
%% Luong Nguyen 10/23/15
%% INPUT: method_name (EGB, QuadTree, NCut, MShift)
%         corresponding parameters
%         EGB: sigma, k, min
%         NCut: num_seg
%         QuadTree: num_seg
%         MShift: 'SpatialBandWidth'(sbw), 'RangeBandWidth'(rbw), 'MinimumRegionArea'(mra)
function params = extract_seism_par(method_name)
seism_dir = '/Users/lun5/Research/github/seism'; %mac
addpath(genpath(seism_dir));
egb_par_dir = fullfile(seism_dir, 'results',method_name);
file_list = dir(fullfile(egb_par_dir,'test_bce*'));
file_list = {file_list.name}';
num_par_settings = length(file_list);
switch method_name
    case 'EGB' 
        num_par = 3;
    case 'MShift'
        num_par = 3;
    case 'NCut'
        num_par = 1;
    case 'QuadTree'
        num_par = 1;
    otherwise
        error(['Unexpected method name. Only input one of the following:' ...
        'EGB, NCut, QuadTree, MShift']);
end
params = zeros(num_par_settings,num_par);

for i = 1:num_par_settings
    [~,fname,~] = fileparts(file_list{i});
    fname_split = strsplit(fname,'_');
    for j = 1:num_par
        params(i,j) = str2double(fname_split{2+j});
    end
end

save(['params_seism' '_' method_name '.mat'],'params')