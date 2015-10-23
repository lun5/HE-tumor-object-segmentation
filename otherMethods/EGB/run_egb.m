%% script to run Felz-Hutt EGB method
% Luong Nguyen 09/21/2015

github_dir = '/home/lun5/github/HE-tumor-object-segmentation';
addpath(genpath(github_dir));
egb_dir = '/home/lun5/Documents/segment';
cd(egb_dir);

im_dir = '/home/lun5/HEproject/data/Tiles_512_ppm';
result_dir = '/home/lun5/HEproject/evaluation_results/EGB';
seg_result_dir = fullfile(result_dir,'segmentations_seism');

if ~exist(result_dir,'dir')
    mkdir(result_dir)
end

if ~exist(seg_result_dir,'dir')
    mkdir(seg_result_dir);
end

im_list = dir(fullfile(im_dir,'*.ppm'));
im_list = {im_list.name}';
% sigma, k, min
tmp = load('params_seism.mat'); params = tmp.params;
%k_vec = 5050:50:10000;
%sig = 0.5; min_val = 20;
run_times = cell(length(im_list),2);
parfor i = 1:length(im_list)
    im_name = im_list{i}(1:end-4);
    T = tic;
    for j = 1:size(params,1)%length(k_vec)
        sig = params(j,1); k = params(j,2); min_val = params(j,3);
        %k = k_vec(j);
        cmm = ['./segment ' num2str(sig) ' ' num2str(k) ' ' ...
            num2str(min_val) ' ' fullfile(im_dir,[im_name '.ppm']) ' ' ...
            fullfile(seg_result_dir,[im_name '_' num2str(k) '.ppm'])];
        s = evalc_parfor(cmm);
    end
    t = toc(T);
    run_times(i,:)= {im_name,t};
    fprintf('Done with image %s in %.2f s\n',im_name, t);
end

save(fullfile(seg_result_dir,'runtimes.mat'),'run_times');
%% convert from ppm image to segmentation results
matfile_result_dir = fullfile(result_dir,'segmented_images_seism');
if ~exist(matfile_result_dir,'dir')
    mkdir(matfile_result_dir);
end

output_list = dir(fullfile(seg_result_dir,'*.ppm'));
output_list = {output_list.name}';
split_output_list = cellfun(@(s) regexp(s,'[_.]','split'), output_list, 'UniformOutput',false);
split_output_list = cat(1,split_output_list{:});
[C,ia,ic] = unique(split_output_list(:,1));

num_segs_per_im = floor(length(ic)/length(C));

parfor i = 1:length(C)
    T = tic;
    im_name = C{i};
    segs = cell(num_segs_per_im,1);
    outFile = fullfile(matfile_result_dir,[im_name,'.mat']);
    if exist(outFile,'file')
        fprintf('Already calculated for image %s\n',imname);
        continue;
    end
    for j = ((i-1)*num_segs_per_im + 1) : i*num_segs_per_im
        check_im_name = strcmp(im_name,split_output_list{j});
        if ~check_im_name
            error('Imnames do not match');
        end
        I = imread(fullfile(seg_result_dir,output_list{j}));
        seg = rgb2label(I);
        segs{j-(i-1)*num_segs_per_im} = seg;
    end
    parsave(outFile,segs);
    t = toc(T);
    fprintf('Done with image %s in %.2f seconds\n',im_name,t);
end

disp('Done');

