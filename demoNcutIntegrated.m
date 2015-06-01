% Demo file
% -------------------------------------------------------------------------
% HE tumor object segmentation 
% Luong Nguyen, lun5@pitt.edu
% Please email me if you find bugs, or have suggestions or questions

function demoNcutIntegrated

%% compile and check for error
%addpath(genpath(pwd));
%tiles_dir = fullfile('Z:\','TilesForLabeling_tiff_renamed'); %window
tiles_dir = '/Users/lun5/Box Sync/TilesForLabeling_tiff_renamed'; %mac
%tiles_dir =  '/home/lun5/HEproject/TilesForLabeling_tiff_renamed'; %linux
clear raw_image Pts ans im mdist opts_affinity opts_clustering which_affinity which_features
%imname = 'uaZFwoHref.tif';
%imname = 'aNaggwovpxANWq0.tif';
%imname = 'jRh62FQ8hUZWlA.tif';
%imname = '0ANZqyIBfUc.tif';
%imname = '95f7k8loesyevi.tif';
%imname = 'cxwrYBYWredN.tif';
%imname = '9uixINHtjjiS.tif';
%imname = 'w8kwtop6hyp.tif';
%imname = '2ALe5NgRyfnpo.tif';
%imname = 'jbaKL4TsEqT.tif';
%imname = 'k2yxq1tbr6kpny0.tif';
%imname = 'vmp8mdxkod3xxzu.tif';
%imname = 'h1402uhfkz.tif';
imname = 'dRfMkOErZY.tif';
%imname = 'ycivjoxn14stvq.tif';
%imname = 'fFwTGXYlhYNa.tif';
%imname = 'pLYZEV43nHWmUDK.tif';
%imname = 'LLV232_D04_20x_max_proj.tif';
%tiles_dir = fullfile(pwd,'test_images');
%imname = '253027.jpg';
%% result directory
% splitStr = regexp(imname,'\.','split');
% imresult_dir = fullfile(pwd,'results','HE_results',[splitStr{1} 'crop2']);
% 
% if ~exist(imresult_dir,'dir')
%     mkdir(imresult_dir);
%     fileattrib(imresult_dir,'+w');
% end
imname = lower(imname);
raw_image = imread(fullfile(tiles_dir, imname));
ndown = 4;dz_im = raw_image(1:ndown:end,1:ndown:end,:);
%figure; imshow(raw_image);rect = getrect;dz_im = imcrop(raw_image,rect);
I = double(dz_im);figure; imshow(I/255);size(I)

opts = setEnvironment_affinity;
which_features = {'hue opp'};
hue_image = getFeatures(I,1,which_features,opts);
hue_image = hue_image{1};
nbSegments = 100; % number of eigen vectors
[SegLabel,NcutDiscrete,NcutEigenvectors,NcutEigenvalues,W,imageEdges]= NcutImage(hue_image,nbSegments);

%% Segment image
% builds an Ultrametric Contour Map from the detected boundaries (E_oriented)
% then segments image based on this map    
% this part of the code is only supported on Mac and Linux    
if (~ispc) && opts.calculate_segments  
    tic;thresh = 0.2;
    E_ucm = contours2ucm_crisp_boundaries(mat2gray(E_oriented));
    segmented_image = ucm2colorsegs(E_ucm,I,thresh);toc;
    if opts.plot_results, 
        figure; subplot(121); imshow(uint8(I)); subplot(122); 
        imshow(uint8(segmented_image)); 
    end;
else
    segmented_image = [];
end

end
