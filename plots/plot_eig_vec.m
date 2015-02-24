%% plot eigen vectors

figure;
ha = tight_subplot(3,3,[.01 .0],[0 0],[0 0]);
for i =1:size(E_oriented,3)
    axes(ha(i));imagesc(E_oriented(:,:,i));
    axis equal; axis tight; axis off; %colormap('gray');
end

figure; imagesc(E);
axis equal; axis tight; axis off; colormap('gray');
