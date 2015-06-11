%% make final image montage
DATA_DIR = '/Users/lun5/Research/data';
IMG_DIR = '/Users/lun5/Box Sync/TilesForLabeling_tiff_renamed';
gt_dir = fullfile(DATA_DIR,'gt_1024_diplay');
orig_img_dir = fullfile(DATA_DIR,'eval_hue_512_512_sig3_exp2_newAffinity_Merge/image1024');
luminance_dir = fullfile(DATA_DIR,'eval_luminance/segmented_images');
col_var_dir = fullfile(DATA_DIR,'eval_col_var_512_sig3_exp2/segmented_images');
hue_dir = fullfile(DATA_DIR,'eval_hue_512_512_sig3_exp2_newAffinity_Merge/segmented_images');

% img_list = {'2ALe5NgRyfnpo.tif','7vj4ekusieek6ys.tif','hrlxfimwgjas.tif',...
%     'fFwTGXYlhYNa.tif','95f7k8loesyevi.tif','9uixINHtjjiS.tif','0ANZqyIBfUc.tif',...
%     'aNaggwovpxANWq0.tif','jbaKL4TsEqT.tif','k1boslrxx7.tif','ocmmhhrtzz5.tif',...
%     'n2wolhpsak70anw.tif','oasyczs5983.tif','p1edtnwreykh.tif',...
%     'ZcU7vNvnWSrKAf.tif','zhB3vhJ2vc.tif','ZleLRgfLYYfT4C.tif',...
%     'zpD2rLbjZm.tif','zXvCwQJOEyD.tif','s5rkUDO7teXwKl.tif'};

IMG_EXT = '.tif';
img_list = dirrec(IMG_DIR,IMG_EXT);
parfor i = 1:length(img_list)
    [~,im_name,~] = fileparts(img_list{i});
    %im_name = lower(im_name);
    filenames = {fullfile(orig_img_dir,[im_name '_512_512.tif']),...
        fullfile(luminance_dir,[im_name '_segmentedImage.tif']),...
        fullfile(col_var_dir,[im_name '_segmentedImage.tif']),...
        fullfile(hue_dir,[im_name '_segmentedImage.tif']),...
        fullfile(gt_dir,[im_name '.tif'])};
    a = cell(5,1);
    for j = 1:length(filenames)
        a{j} = imread(filenames{j});
    end
    montage = [a{1} a{2} a{3} a{4} a{5}];
    imwrite(montage,fullfile(DATA_DIR,'montageImageSegment',[im_name '_montage.tif']),'Resolution',300);
    %
    %         %filenames = {img_list{i}, fullfile(RESULTS_DIR,'segmented_images',[im_name '_segmentedImage.tif'])};
    %         if ~exist(fullfile(DATA_DIR,'montageImageSegment',[im_name '_montage.tif']),'file')
    %         h = montage(filenames,'Size',[1 5]);
    %         MyMontage = getframe(gca); %//
    %         imwrite(MyMontage.cdata,fullfile(DATA_DIR,'montageImageSegment',[im_name '_montage.tif']),'Resolution',300);
             fprintf('\n Print montage image of %s...',im_name); close all;
    %         end
end