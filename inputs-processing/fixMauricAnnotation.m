% script to fix Maurice's annotation
result_dir = 'Z:\HEproject\evaluation_results\eval_non_expert\Maurice';
imname = 'sp2byg1b33ghl';
im_list = {'1bhjqxecct','3tdnf6bhanan','4nkj5wqcqj',...
    '9fukmjujwi','13nedzdzfh','apaal7fc2payi','aqizfuqbbxyu','bylklqnsvf4d',...
    'drfmkoerzy','finkidqlnznihk','hrdqmlu2ig','ovvp7sjwtrrax3t','sfzywlg4291ch9',...
    'shyy81fuot9hu','sp2byg1b33ghl','t42x5lgm4py','taljyo23jlxd','tj6ac55wqwfi',...
    'uhin9nl4ju7bls','updte6rxh7afsf','uraxeh1spli7ky9', 'uyuukdfzqq','vdqmu8xq2xc',...
    'vm3qo9caekfodi','vwxsdt6g3n2f','vyggwhmvuj1','w8kwtop6hyp','w9fpyfgxtibvdk',...
    'wm9g9oyreu','xhcy7eroyqz','ycivjoxn14stvq','ylf5dghxyp2sxp','ynqnt1mljk4e',...
    'yxr1wm9fonpvmf','yxt9szgtneh8','z8wtvu0v8g','zlelrgflyyft4c','zxvcwqjoeyd'};
if ~exist(fullfile(result_dir,'segmented_images_new'),'dir')
    mkdir(fullfile(result_dir,'segmented_images_new'));
end
for i = 1:length(im_list)
    imname = im_list{i};
    tmp = load(fullfile(result_dir,'segmented_images',[imname '.mat']));
    seg = tmp.data{1};
    figure; imshow(label2rgb(seg));
    seg = seg(1:128,1:128);
    seg = uint8(imresize(seg,4));
    figure; imshow(label2rgb(seg)); 
    parsave(fullfile(result_dir,'segmented_images_new',[imname '.mat']),{seg}),
end