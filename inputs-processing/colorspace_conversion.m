csvwrite(fullfile(pwd,'results','tp10-611gland7.csv'),theta);
csvwrite(fullfile(pwd,'results','gland1_oppCol.csv'),f_maps(:)');
im_rgb = double(I);
im_lab = colorspace('lab<-rgb',im_rgb./255); range(im_lab) 
im_lch = colorspace('lch<-rgb',im_rgb./255); range(im_lch)
hue_lch = im_lch(:,:,3); hue_rad = deg2rad(hue_lch(:));
hue_rad_im = deg2rad(hue_lch) - pi;
imtool(hue_rad_im,[]);
figure; rose(hue_rad);

%csvwrite(fullfile('results','tp10-611_gland1snip_theta_oppCol.csv'),theta);
%csvwrite(fullfile(pwd,'results','tp10-611gland7_theta_lch.csv'),hue_rad');
csvwrite(fullfile(pwd,'results','gland1_lch.csv'),hue_rad');
[h,s,v] = rgb2hsv(im_rgb./255); range(h)
figure; rose(hue_rad); % makes sense since it is between red and blue

csvwrite(fullfile(pwd,'results','gland3_snip_opp.csv'),theta);

im_nature = imread(fullfile(pwd,'test_images','48055.jpg'));
r_nature = im_nature(:,:,1);g_nature = im_nature(:,:,2);b_nature = im_nature(:,:,3);
im_nature_rgb = cat(1,r_nature(:)',g_nature(:)',b_nature(:)');
ind_perm = randperm(numel(r_nature));
[U,~,~] = svd(double(im_nature_rgb(:,ind_perm(1:50000))),0);
rotation_matrix = [-U(:,1) U(:,2:3)]';
