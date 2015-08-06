%% Process output of TextonBoost to evaluate its performance
% Luong Nguyen 08/06/2015

% Read in the files
% input_dir = 'Z:\TextonInput\Temp';
% input_fnames = dir(fullfile(input_dir,'*.png'));
% input_fnames = {input_fnames.name}';
% num_files = length(input_fnames);
% 
% %% change the image names
% parfor i = 1:num_files
%     fname = input_fnames{i};
%     [~,im_name,~] = fileparts(fname);
%     im_name = strsplit(im_name,{'-' '.'});
%     im_name = im_name{end-1};
%     movefile(fullfile(input_dir,fname),fullfile(input_dir,[im_name,'.png']));
% end

%% convert image --> labels
input_dir = 'Z:\TextonInput\Temp';
input_fnames = dir(fullfile(input_dir,'*.png'));
input_fnames = {input_fnames.name}';
num_files = length(input_fnames);
components = {'carcinoma','vessel','fat','stroma','lymphocyte','duct',...
    'white space','atypical ductal hyperplasia'};%,'benign terminal lobular unit','unsure'};
num_comps = length(components);
color_map = zeros(num_comps,3);
% mapping colors
for id = 1: length(components)
   color_map(id,:) = mapLabel2Color(id-1);
end
output_dir = fullfile('Z:\TextonInput\output');
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
parfor i = 1:num_files
    [~,im_name,~] = fileparts(input_fnames{i});
    im = imread(fullfile(input_dir,[im_name,'.png']));
    labels = mapColor2Label(reshape(im,[size(im,1)*size(im,2) size(im,3)]));
    seg = reshape(labels,[size(im,1) size(im,2)]); %segmentation
    bmap = logical(seg2bdry(seg,'imageSize'));  
    parsave(fullfile(output_dir,[im_name,'_seg.mat']),seg);
    parsave(fullfile(output_dir,[im_name,'_bdry.mat']),bmap);
    %figure; imshow(mat2gray(bmap));
    fprintf('Done with file %s\n',im_name);
end

%% evaluate F-score and such
evDir = fullfile(out_dir,'evFiles');
if ~exist(evDir,'dir')
    mkdir(evDir);
end
thinpb = true;
maxDist = 0.01; % 0.0075
nthresh = 99; 
gtDir = 'Z:\HEproject\data\groundTruth_512_512';
parfor i = 1:num_files,
    T = tic;
    [~,im_name,~] = fileparts(input_fnames{i});
    fprintf('\n\nEvaluate results %s...',im_name);
    evFile4 = fullfile(evDir, strcat(im_name, '_ev4.txt'));
    if ~isempty(dir(evFile4)), continue; end
    
    inFile = fullfile(output_dir, strcat(im_name, '_seg.mat'));
    gtFile = fullfile(gtDir, strcat(im_name, '.mat'));
    evFile1 = fullfile(evDir, strcat(im_name,'_ev1.txt'));
    evFile2 = fullfile(evDir, strcat(im_name, '_ev2.txt'));
    evFile3 = fullfile(evDir, strcat(im_name, '_ev3.txt'));

    evaluation_bdry_image_seg(inFile,gtFile, evFile1, nthresh, maxDist, thinpb);
    evaluation_reg_image_seg(inFile, gtFile, evFile2, evFile3, evFile4, nthresh);    
    %disp(i);
    t = toc(T); fprintf('done: %1.2f sec\n', t);
end
collect_eval_bdry(evDir);
collect_eval_reg(evDir);
[ Fop_ods, P_ods, R_ods, bestT, Fop_ois, P_ois, R_ois] = eval_Fop(input_dir, gtDir, output_dir, evDir);