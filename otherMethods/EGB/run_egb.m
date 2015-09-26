%% script to run Felz-Hutt EGB method
% Luong Nguyen 09/21/2015

github_dir = '/home/lun5/github/HE-tumor-object-segmentation';
addpath(genpath(github_dir));
egb_dir = '/home/lun5/Documents/segment';
cd(egb_dir);

im_dir = '/home/lun5/HEproject/data/Tiles_512_ppm';
result_dir = '/home/lun5/HEproject/evaluation_results/EGB';
seg_result_dir = fullfile(result_dir,'segmentations');

if ~exist(result_dir,'dir')
    mkdir(result_dir)
end

if ~exist(seg_result_dir,'dir')
    mkdir(seg_result_dir);
end

im_list = dir(fullfile(im_dir,'*.ppm'));
im_list = {im_list.name}';
k_vec = 100:50:2000;
sig = 0.5; min_val = 20;
run_times = cell(length(im_list),2);
parfor i = 1:length(im_list)
    im_name = im_list{i}(1:end-4);
    T = tic;
    for j = 1:length(k_vec)
        k = k_vec(j);
        cmm = ['./segment ' num2str(sig) ' ' num2str(k) ' ' ...
            num2str(min_val) ' ' fullfile(im_dir,[im_name '.ppm']) ' ' ...
            fullfile(seg_result_dir,[im_name '_' num2str(k) '.ppm'])];
        s = evalc_parfor(cmm);
    end
    t = toc(T);
    run_times(i,:)= {im_name,t};
    fprintf('Done with image %s in %.2f s\n',im_name, t);
end

save(fullfile(result_dir,'runtimes.mat',run_times));
%% convert from ppm image to segmentation results
matfile_result_dir = fullfile(result_dir,'matfiles');
if ~exist(matfile_result_dir,'dir')
    mkdir(matfile_result_dir);
end

output_list = dir(fullfile(seg_result_dir,'*.ppm'));
output_list = {output_list.name}';
split_output_list = cellfun(@(s) regexp(s,'[_.]','split'), output_list, 'UniformOutput',false);
split_output_list = cat(1,split_output_list{:});
[C,ia,ic] = unique(split_output_list(:,1));

% I need a dictionary here
%keySet = C;
%valueSet = repmat({{2}},size(C));
seg_container = cell(length(split_output_list),1); %containers.Map(keySet,valueSet);
parfor i = 1:length(split_output_list);
    im_name = split_output_list{i,1};
    I = imread(fullfile(seg_result_dir,output_list{i}));
    seg = rgb2label(I);
    seg_container{i} = seg;
    %seg_container(im_name) = cat(1, seg_container(im_name),{seg});    
end

for i =1:length(C)
    segs = seg_container(ic == i);
    parsave(fullfile(matfile_result_dir,[key '.mat']),segs);
end

parfor key = seg_container.keys
    segs = seg_container(key);
    parsave(fullfile(matfile_result_dir,[key '.mat']),segs{2:end});
end

disp('Done');

