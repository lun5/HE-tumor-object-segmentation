im_lab = colorspace('lab<-rgb',im_rgb./255); range(im_lab) 
im_lch = colorspace('lch<-rgb',im_rgb./255); range(im_lch)
hue_lch = im_lch(:,:,3); hue_rad = deg2rad(hue_lch(:));
csvwrite(fullfile(pwd,'results','gland3_snip_lch.csv'),hue_rad');
[h,s,v] = rgb2hsv(im_rgb./255); range(h)
figure; rose(hue_rad); % makes sense since it is between red and blue

csvwrite(fullfile(pwd,'results','gland3_snip_opp.csv'),theta);