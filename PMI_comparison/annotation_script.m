GITHUB = 'C:\Users\luong_nguyen\Documents\GitHub';
HESEG = fullfile(GITHUB,'HE-tumor-object-segmentation');
cd(HESEG); addpath(genpath(HESEG));
LABELME = fullfile(GITHUB,'LabelMeToolBox');
addpath(genpath(LABELME)); % we could put this into the toolbox of H&E
HOMEIMAGES = fullfile(LABELME,'Images'); % you can set here your default folder
HOMEANNOTATIONS = fullfile(LABELME,'Annotations'); % you can set here your default folder
% select one annotation file from one of the folders:
filename = fullfile(HOMEANNOTATIONS, 'renamed_images', '4d0ylpdlwf.xml');
% This line reads the entire database into a Matlab struct
%database = LMdatabase(fullfile(HOMEANNOTATIONS,'renamed_images'));
% read the image and annotation struct:
[annotation, img] = LMread(filename, HOMEIMAGES);
folderlist = {'renamed_images'};%filelist = {'9uixinhtjjis.jpg'};
filelist = {'5aoqp0sbfxswmz.jpg'};
database = LMdatabase(HOMEANNOTATIONS, HOMEIMAGES, folderlist, filelist);
[annotation, img] = LMread(database,1, HOMEIMAGES);
% plot the annotations
LMplot(annotation, img)

[D_tumor,j_tumor] = LMquery(database, 'object.name', 'tumor');
[D_stroma,j_stroma] = LMquery(database, 'object.name', 'stroma');
[D_lymp,j_lymp] = LMquery(database, 'object.name', 'chronic inflammation');
LMdbshowscenes(D_tumor, HOMEIMAGES);
LMdbshowscenes(D_stroma, HOMEIMAGES);
%% tumor
dd = D_tumor(1);
dd = D_stroma(1);
%dd = D_lymp(1);
Nobjects = length(dd.annotation.object);n =0;

% for i = 1:length(D_tumor)
%     LMdbshowobjects(D_tumor(i),HOMEIMAGES);
% end
sourcedir = 'Z:\';
tiles_dir = fullfile(sourcedir,'TilesForLabeling_tiff_renamed');

%img = imread(fullfile(tiles_dir,'9uixinhtjjis.tif')); % Load image
img =imread(fullfile(tiles_dir,'5aoqp0sbfxswmz.tif'));
[nrows ncols c] = size(img);
for j = 1:Nobjects
    n = n+1;
    [X,Y] = getLMpolygon(dd.annotation.object(j).polygon);
    
    BW = poly2mask(double(X),double(Y), nrows, ncols);
    polygonRegion = img.*repmat(uint8(BW),[1 1 3]);
    
    crop(1) = max(min(X)-2,1);crop(2) = min(max(X)+2,ncols);
    crop(3) = max(min(Y)-2,1);crop(4) = min(max(Y)+2,nrows);
    crop = round(crop);
                
    % Image crop:
    imgCrop = polygonRegion(crop(3):crop(4), crop(1):crop(2), :);

    figure; imshow(imgCrop)

BW_crop = BW(crop(3):crop(4),crop(1):crop(2));
I = double(imgCrop);
opts = setEnvironment_affinity;
Nsamples = 10000;
which_features = opts.features.which_features;
f_maps = getFeatures(I,1,which_features,opts);
[F,p1,p2] = sampleF(f_maps{1},Nsamples,opts,BW_crop);

figure; imshow(uint8(BW_crop)*255); hold on;
plot(p1(:,2), p2(:,1),'bx','MarkerSize',10);
plot(p2(:,2), p2(:,1),'rx','MarkerSize',10);
hold off;

Fsym = [F; F(:,2) F(:,1)];
figure; ndhist(Fsym(:,1),Fsym(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
xlabel('\phi'); ylabel('\psi'); axis square;set(gca,'FontSize',16);
set(gcf,'color','white') 

end
% [ params,~, prior_probs] = mixture_of_bivariate_VM(F, 6);
% mixture_params.params = params;
% mixture_params.prior_probs = prior_probs;
% plotPMI_theta;
