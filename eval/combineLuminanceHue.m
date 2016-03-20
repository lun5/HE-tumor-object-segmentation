% generate E_oriented files in the last desperate attempt
% Luong Nguyen
% 3/16/16

hue_dir = '/Users/lun5/Research/data/normalized_evaluation_results/PMI_lowres_accurate/E_oriented/';
luminance_dir = '/Users/lun5/Research/data/normalized_evaluation_results/Isola_speedy/E_oriented/';

imlist = dir(fullfile(hue_dir,'*.mat'));
imlist = {imlist.name}';

out_dir = '/Users/lun5/Research/data/normalized_evaluation_results/combined_1_1';
if ~exist(out_dir,'dir');
    mkdir(out_dir);
    mkdir(fullfile(out_dir,'E_oriented'));
    mkdir(fullfile(out_dir,'edge_map'));
end

parfor i = 1:length(imlist)
    imname = imlist{i}(1:end-15);
    fprintf('start with image %s ...',imname); T = tic;
    tmp = load(fullfile(hue_dir,[imname '_E_oriented.mat']));
    E_oriented_hue = tmp.data;
    tmp = load(fullfile(luminance_dir,[imname '_E_oriented.mat']));
    E_oriented_luminance = tmp.data;
    E_oriented = E_oriented_hue + E_oriented_luminance;
    %E_oriented = (E_oriented - min(E_oriented(:)))./(max(E_oriented(:)) - min(E_oriented(:)));
    E = max(E_oriented,[],3);
    imwrite(mat2gray(1-E),fullfile(out_dir,'edge_map',[imname '.tif']));
    parsave(fullfile(out_dir,'E_oriented',[imname '_E_oriented.mat']),E_oriented);
    fprintf('Done in %d seconds\n',toc(T));
end
    
