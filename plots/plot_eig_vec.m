%% plot eigen vectors

figure;
for i =1:size(E_oriented,3)
    subplot(4,2,i); imagesc(E_oriented(:,:,i));
    axis equal;axis([0 im_sizes{1}(1) 0 im_sizes{1}(2)])
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    %colormap('gray');
end

figure; imagesc(E);
axis equal;axis([0 im_sizes{1}(1) 0 im_sizes{1}(2)])
set(gca,'xtick',[]);set(gca,'ytick',[]);
colormap('gray');
