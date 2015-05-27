set(gcf,'PaperPositionMode','auto')
print(fullfile('results','pmi'),'-dtiff','-r300');

set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')
print(fullfile('results','Ws_loc'),'-dtiff','-r300');

set(gcf,'color','white');
colorbar off;
set(gcf,'PaperPositionMode','auto')
print(fullfile('results','Ws_1'),'-dtiff','-r300');
set(gcf,'PaperPositionMode','auto')
print(fullfile('results','Ws_2'),'-dtiff','-r300');
set(gcf,'PaperPositionMode','auto')
print(fullfile('results','Ws_3'),'-dtiff','-r300');

imwrite(I/255,fullfile('results','im.tif'),'Resolution',300);
imwrite(1-mat2gray(E),fullfile('results','edges.tif'),'Resolution',300);
imwrite(uint8(segmented_image),fullfile('results','segmented.tif'),'Resolution',300);