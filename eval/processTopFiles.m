DATA_DIR ='/home/lun5/HEproject/'; % linux
IMG_DIR = '/Users/lun5/Box Sync/TilesForLabeling_tiff_renamed'; %mac
GT_DIR = '/Users/lun5/Research/data/groundTruth_512_512';
RESULTS_DIR = '/Users/lun5/Research/data/eval_hue_512_512_sig3';
outDir = fullfile(RESULTS_DIR,'outputTop');
if ~exist(outDir,'dir')
    mkdir(outDir)
end

eval_bdry_img = dlmread(fullfile(RESULTS_DIR,'eval_bdry_img.txt'));
ind_top = find(eval_bdry_img(:,5) > 0.4);
IMG_EXT = '.tif';
img_list = dirrec(IMG_DIR,IMG_EXT);
img_list = img_list(:,ind_top);

thinpb = true; 
maxDist = 0.01; 
nthresh = 99;

parfor i = 1:numel(img_list),
    tic;
    [~,im_name,~] = fileparts(img_list{i});
    evFile4 = fullfile(outDir, strcat(im_name, '_ev4.txt'));
    if ~isempty(dir(evFile4)), continue; end
    
    inFile = fullfile(RESULTS_DIR, strcat(im_name, '.mat'));
    gtFile = fullfile(GT_DIR, strcat(im_name, '.mat'));
    evFile1 = fullfile(outDir, strcat(im_name,'_ev1.txt'));
    evFile2 = fullfile(outDir, strcat(im_name, '_ev2.txt'));
    evFile3 = fullfile(outDir, strcat(im_name, '_ev3.txt'));

    evaluation_bdry_image(inFile,gtFile, evFile1, nthresh, maxDist, thinpb);
    evaluation_reg_image(inFile, gtFile, evFile2, evFile3, evFile4, nthresh);    
    %disp(i);
    fprintf('\n\nEvaluate results %s...',im_name);
    t = toc; fprintf('done: %1.2f sec\n', t);
end

collect_eval_bdry(outDir);
collect_eval_reg(outDir);
eval_Fop(GT_DIR, RESULTS_DIR, outDir);
plot_eval(outDir);

