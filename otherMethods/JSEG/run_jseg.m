%% Run JSEG algorithm on the data
% This actually behaves quite well
% Luong Nguyen 10/05/2015

github_dir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation';
jseg_dir = fullfile(github_dir,'otherMethods','JSEG');
cd(jseg_dir);
im_dir = 'Z:\Tiles_512_jpg';
%output_dir = fullfile('Z:\HEproject\evaluation_results\JSEG','multi_scale');
output_dir = fullfile('Z:\HEproject\evaluation_results\JSEG','one_scale');

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

if ~exist(fullfile(output_dir,'seg_im'),'dir')
    mkdir(fullfile(output_dir,'seg_im'));
end

if ~exist(fullfile(output_dir,'bdry_im'),'dir')
    mkdir(fullfile(output_dir,'bdry_im'));
end

if ~exist(fullfile(output_dir,'mat_files'),'dir')
    mkdir(fullfile(output_dir,'mat_files'));
end

if ~exist(fullfile(output_dir,'gif_files'),'dir')
    mkdir(fullfile(output_dir,'gif_files'));
end

quantize_threshold = 0:25:600;
num_thres = length(quantize_threshold);
im_list = dir(fullfile(im_dir,'*.jpg'));
im_list = {im_list.name}';
num_images = length(im_list);
quote = '''';
for i = 1:num_images
    T = tic; 
    im_name = im_list{i}(1:end-4);
    fprintf('Start with image %s...',im_name);
    for j = 1: num_thres
        q_thresh = quantize_threshold(j);
        gif_file = fullfile(output_dir,'gif_files',[im_name '_'  num2str(q_thresh) '.gif']);
        if ~exist(gif_file,'file')
            %expr = ['segwin -i ', fullfile(im_dir,im_list{i}), ' -t 6 -o ', ...
            %  fullfile(output_dir,[im_name '.jpg']), ' 0.9 -r9 ', ...
            %  fullfile(output_dir,[im_name '.gif']), ' -q 255'];
            expr = ['segwin -i ', fullfile(im_dir,im_list{i}), ' -t 6 -r9 ', ...
                gif_file,' -l 1 -q ' num2str(q_thresh)];
            out_expr = evalc(['system(' quote expr quote ')']);
        end
    end
    t = toc(T); fprintf(' Done in %.2f seconds\n',t);
end

%quantize_threshold = 50:50:600;
%num_thres = length(quantize_threshold);

parfor i = 1:num_images
    T = tic; 
    im_name = im_list{i}(1:end-4);
    fprintf('Start with image %s...',im_name);
    I = imread( fullfile(im_dir,im_list{i}));
    segs = cell(num_thres,1);
    for j = 1:num_thres
        q_thresh = quantize_threshold(j);
        bdry_im_fname = fullfile(output_dir,'bdry_im',[im_name, '_' num2str(q_thresh), '_bdry.jpg']);
        labels = imread(fullfile(output_dir,'gif_files',[im_name '_'  num2str(q_thresh) '.gif']),1);
        segs{j} = labels+1;
        if ~exist(bdry_im_fname,'file')
            edge_map = seg2bdry(labels,'imageSize');
            % change thickness of edges
            edge_map = imdilate(edge_map, strel('disk',1));
            edge_map_im = I.*uint8(repmat(~edge_map,[1 1 3]));
            imwrite(edge_map_im,fullfile(output_dir,'bdry_im',[im_name, '_' num2str(q_thresh), '_bdry.jpg']));
            imwrite(label2rgb(labels),fullfile(output_dir,'seg_im',[im_name, '_' num2str(q_thresh), '_seg.jpg']));
        end
    end
    parsave(fullfile(output_dir,'mat_files',[im_name '.mat']),segs);
    t = toc(T); fprintf(' Done in %.2f seconds\n',t);
end

disp('done');
%color quantization threshold - specify  
%values 0-600, leave blank for automatic determination. 
%The higher the value, the less number of quantized colors in the image. 
%For color images, try 250. If you are unsatisfied with the result 
%because two neighboring regions with similar colors are not getting separated, 
%please try a smaller value say 150.