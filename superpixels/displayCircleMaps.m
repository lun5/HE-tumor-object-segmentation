%% Turn the circle into display
%% Luong Nguyen 04/16/2015
clearvars;
sourcedir = 'Z:\';
superpixel_input_dir = fullfile(sourcedir,'super_pixels');

circleFnames = dir(fullfile(superpixel_input_dir,'*_se1_minNuc3_minStr5_minLum10_circle_map'));
imagepaths = {circleFnames.name}';
imagepaths = {'036YaafqbAvH5su_se1_minNuc3_minStr5_minLum10_circle_map',...
    '0ANZqyIBfUc_se1_minNuc3_minStr5_minLum10_circle_map'};
numImages = length(imagepaths);% 420
%cd(superpixel_input_dir)
M=[.41 .10 .310 ;  1 0.31 0.8 ; 0.31 0.9 1 ; 0 1 1 ;...
    1 1 0 ; 1 0 1 ; 1 1 1 ; 0.5 0 0 ; 0 0 0.5 ];

for j = 1:numImages
    ObjectsFileName = imagepaths{j}; 
    d1 = dlmread(ObjectsFileName, ',', 1, 1);
    b=label2rgb(d1,M);
    im_splitStr = regexp(ObjectsFileName,'\_','split');
    fname = fullfile(superpixel_input_dir,[im_splitStr{1},'_se1_minNuc3_minStr5_minLum10_display_circle.jpg']);
%     figure, imshow(b);
%     set(gca,'LooseInset',get(gca,'TightInset'))
%     set(gcf,'color','white') % White background for the figure.
    imwrite(b,fname);
    sprintf('Finish with file %s\n',im_splitStr{1})
end
