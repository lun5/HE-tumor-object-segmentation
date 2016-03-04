DATA_DIR = 'Z:\';
IMG_DIR = 'Z:\Tiles_512\Test';
GT_DIR = 'Z:\HEproject\data\groundTruth_512_512';
RESULTS_DIR = fullfile(DATA_DIR,'HEproject','evaluation_results','eval_non_expert','Maurice');
RESULTS_DIR = fullfile(DATA_DIR,'HEproject','evaluation_results','eval_non_expert','Om');
SEG_DIR = fullfile(RESULTS_DIR,'segmented_images');
OUT_DIR = fullfile(RESULTS_DIR,'bdry_disp');
if ~exist(OUT_DIR,'dir');
    mkdir(OUT_DIR);
end
img_list = dirrec(SEG_DIR,'.mat');
num_im = length(img_list);
parfor i = 1:num_im
    [~,im_name,~] = fileparts(img_list{i}); im_name = lower(im_name);
    tmp = load(fullfile(SEG_DIR,[im_name '.mat']));
    I = imread(fullfile(IMG_DIR,[im_name '.tif']));
    outFile = fullfile(OUT_DIR,[im_name, '.tif']);
    segs = tmp.data;
    edge_map = edge(segs{1});
    edge_map = imdilate(edge_map, strel('disk',1));
    edge_map_im = I.*uint8(repmat(~edge_map,[1 1 3]));
    imwrite(edge_map_im,outFile,'Resolution',300);
end
