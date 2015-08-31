% This script is to evaluate the smoothness of signals across H&E images
% the smoothness is increased using a 2D disk filter
% Luong Nguyen 8/18/15

indx = 250;
% plot rgb
figure; plot(im_rgb(indx,:,1),'r-'); hold on;
plot(im_rgb(indx,:,2),'g-');
plot(im_rgb(indx,:,3),'b-');axis tight
ylim([0 1]); set(gca,'FontSize',20);

% plot hsv
hsv = rgb2hsv(X);
hsv = reshape(hsv,size(im_rgb));
figure; plot(hsv(indx,:,1)); hold on
%plot(hsv(indx,:,2));plot(hsv(indx,:,3));
axis tight; ylim([0 1]); set(gca,'FontSize',20);
% hue H&E
figure; plot(im(indx,:)); 
axis tight; ylim([-pi pi]);set(gca,'FontSize',20);

% plot lab
lab = rgb2lab(im_rgb);
figure; plot(lab(indx,:,1)); hold on
plot(lab(indx,:,2));
plot(lab(indx,:,3));axis tight;
set(gca,'FontSize',20);
legend('L','a','b');
% 3 channels
im_1 = reshape(rotated_coordinates(1,:),[size(im_rgb,1) size(im_rgb,2)]);
im_2 = reshape(rotated_coordinates(2,:),[size(im_rgb,1) size(im_rgb,2)]);
im_3 = reshape(rotated_coordinates(3,:),[size(im_rgb,1) size(im_rgb,2)]);
figure; plot(im_1(indx,:),'r-');hold on;axis tight;
plot(im_2(indx,:),'g-');plot(im_3(indx,:),'b-');
set(gca,'FontSize',20);
legend('br','sat*cos(h)', 'sat*sin(h)');

% filter
h = fspecial('disk',5);
I = imfilter(raw_image,h);
% straight line
figure; imshow(im_rgb);
axis on; set(gca,'FontSize',20);
hold on; plot(1:size(im_rgb,1), ones(size(im_rgb,1),1)*indx,'k-','LineWidth',3)