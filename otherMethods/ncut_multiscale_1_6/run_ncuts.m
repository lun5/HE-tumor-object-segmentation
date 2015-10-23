%% script to run NCut method
% Luong Nguyen 09/21/2015
% Ncut multiscale code

% github_dir = '/home/lun5/github/HE-tumor-object-segmentation';
% ncut_dir = fullfile(github_dir,'otherMethods','ncut_multiscale_1_6');
% im_dir = '/home/lun5/HEproject/data/Tiles_512';
% result_dir = '/home/lun5/HEproject/evaluation_results/ncut_multiscale_1_6';

githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; % window
ncut_dir = fullfile(github_dir,'otherMethods','ncut_multiscale_1_6');
im_dir = 'Z:\Tiles_512';
result_dir = fullfile('Z:\HEproject','evaluation_results','ncut_multiscale_1_6');

% github_dir = '/Users/lun5/Research/github/HE-tumor-object-segmentation'; % mac
% ncut_dir = fullfile(github_dir,'otherMethods','ncut_multiscale_1_6');
% im_dir = '/Users/lun5/Research/data/Tiles_512';
% result_dir = '/Users/lun5/Research/data/evaluation_results/ncut_multiscale_1_6';

addpath(genpath(github_dir));
addpath(genpath(ncut_dir));

matfile_result_dir = fullfile(result_dir,'segmented_images');

if ~exist(result_dir,'dir')
    mkdir(result_dir)
end

if ~exist(matfile_result_dir,'dir')
    mkdir(matfile_result_dir);
end

if ~exist(fullfile(result_dir,'seg_im'),'dir')
    mkdir(fullfile(result_dir,'seg_im'));
end

if ~exist(fullfile(result_dir,'bdry_im'),'dir')
    mkdir(fullfile(result_dir,'bdry_im'));
end

im_list = dir(fullfile(im_dir,'*.tif'));
im_list = {im_list.name}';
% sigma, k, min
tmp = load(fullfile(github_dir,'otherMethods','params_seism_NCut.mat')); params = tmp.params;
max_nSegs = max(params);
num_segs = size(params,1);
%run_times = cell(length(im_list),2);
parfor i = 1:8%length(im_list)
    im_name = im_list{i}(1:end-4);
    I = imread(fullfile(im_dir,im_list{i})); 
    segs = cell(num_segs,1);
    outFile = fullfile(matfile_result_dir,[im_name,'.mat']);
    if ~exist(outFile,'file')
        T = tic;
        % calculate the eigenvectors
        outFile_max_nSegs = fullfile(matfile_result_dir,[im_name,'_maxSegs.mat']);
        if ~exist(outFile_max_nSegs,'file')
            [p,q,r] = size(I);
            [layers,C]=compute_layers_C_multiscale(p,q);
            dataW = computeParametersW(I);
            %t =cputime;
            W=computeMultiscaleW(double(I),layers,dataW,[]);
            %fprintf('Calculate %d eigenvectors in %.2f seconds\n',nSegments,cputime-t);
            if ~isempty(C)
                [X,~,~] = computeNcutConstraint_projection(W,C,max_nSegs);
            else
                [X,~,~] = computeKFirstEigenvectors(W,max_nSegs);
            end
            indPixels = (1:p*q)';
            X = reshape(X(indPixels,:),p,q,max_nSegs);
            parsave(outFile_max_nSegs,X);
        else
            tmp = load(outFile_max_nSegs); X = tmp.data;
        end
        % calculate segments
        for j = 1:size(params,1)
            nseg = params(j);
            [SegLabel,~] =discretisation(X(:,:,1:nseg));
            segs{j} = SegLabel;
            bdry_fname = fullfile(result_dir,'bdry_im',[im_name '_' ...
                num2str(nseg) '.jpg']);
            seg_fname = fullfile(result_dir,'seg_im',[im_name '_sbw' ...
                num2str(nseg) '.jpg']);
            if ~exist(bdry_fname,'file')
                edge_map = seg2bdry(segs{j},'imageSize');
                edge_map = imdilate(edge_map, strel('disk',1));
                edge_map_im = I.*uint8(repmat(~edge_map,[1 1 3]));
                imwrite(edge_map_im,bdry_fname);
                imwrite(label2rgb(segs{j}),seg_fname);
            end
        end
           
        parsave(outFile,segs);
        t = toc(T);
        %run_times(i,:)= {im_name,t};
        fprintf('Done with image %s in %.2f s\n',im_name, t);
    end
end

disp('Done');

