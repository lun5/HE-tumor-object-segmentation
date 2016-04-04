%% get the plots for ground truth and non expert observers
DATA_DIR = 'Z:\HEproject';
IMG_DIR = 'Z:\Tiles_512';%\Test';
GT_DIR = fullfile(DATA_DIR,'data','GroundTruth','coarse_fine_GT_512_512');
non_expert_1_dir = 'Z:\HEproject\evaluation_results\eval_non_expert\Maurice';
non_expert_2_dir = 'Z:\HEproject\evaluation_results\eval_non_expert\Om';

input_dirs = {fullfile(GT_DIR,'all_files'),fullfile(non_expert_1_dir,'segmented_images_new'),...
    fullfile(non_expert_2_dir,'segmented_images')};
method_names = {'ground-truth','non-expert','non-expert'};
output_dirs = {fullfile(GT_DIR,'pad_300_bigger'),fullfile(non_expert_1_dir,'pad_300'),...
    fullfile(non_expert_2_dir,'pad_300')};
for i = 1:length(output_dirs)
    if ~exist(output_dirs{i},'dir')
        mkdir(output_dirs{i})
    end
end

for med = 1%2:length(input_dirs)
    img_list = dirrec(input_dirs{med},'.mat');
    for j = 1:length(img_list)
        [~,im_name,~] = fileparts(img_list{j}); im_name = lower(im_name); 
        bdry_outFile = fullfile(output_dirs{med},[im_name, '.tif']);
        I = imread(fullfile(IMG_DIR,[im_name '.tif']));
        tmp = load(fullfile(input_dirs{med},[im_name '.mat']));
        segs = tmp.groundTruth{1}.Segmentation;
        bdry_edge_map = logical(seg2bdry(segs,'imageSize'));
        %segs = tmp.data;
        %bdry_edge_map = logical(seg2bdry(segs{1},'imageSize'));
        bdry_edge_map = imdilate(bdry_edge_map, strel('disk',2));
        bdry_edge_map_im = I.*uint8(repmat(~bdry_edge_map,[1 1 3]));
        %pad_im = padarray(bdry_edge_map_im,[60 60],255,'both');
        %pad_im = insertText(pad_im,[150 5],method_names{med},'FontSize',50,'BoxColor','white');
        pad_im = padarray(bdry_edge_map_im,[300 60],255,'both');
        pad_im = insertText(pad_im,[150 250],method_names{med},'FontSize',50,'BoxColor','white');
        pad_im = pad_im(250:end,:,:);
        imwrite(pad_im,bdry_outFile,'Resolution',300);
    end
end


