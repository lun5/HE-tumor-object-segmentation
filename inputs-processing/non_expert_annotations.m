%% convert non expert annotations to polygon for further processing
IMG_DIR = 'Z:\TilesForLabeling_tiff_renamed';

annot_dir = 'Z:\HEproject\data\GroundTruth\non_expert';
user_name = 'Maurice';
user_name = 'Om';

user_annot_dir = fullfile(annot_dir,user_name);

fname = fullfile(user_annot_dir,'annotations.json');
text = fileread(fname);text = strrep(text,'"','');
%coords_files = strsplit(text,'polygons'); %Maurice
coords_files = strsplit(text,'filename:');
coords_files = coords_files(2:end);
num_files = length(coords_files);
rows = 2048; cols = 2048;
mult = 4;
parfor i = 1:num_files
    T = tic;
    coordinates = regexp(coords_files{i},'coords:','split');
%     fname = regexp(coordinates{end},'filename:','split'); %Maurice
%     coordinates{end} = fname{1};
%     fname = strtrim(regexprep(fname{end},'[},{]',''));
    fname = strsplit(coordinates{1},',');
    fname = strtrim(fname{1});
    coordinates = coordinates(2:end);
    seg =zeros(rows,cols);
    for j = 1:length(coordinates)
        txt = coordinates{j};
        txt = regexprep(txt,'[^a-zA-Z0-9]',' ');
        txt = textscan(txt,'%d');
        coords = double(reshape(txt{1},[2 length(txt{1})/2]));
        mask = poly2mask(coords(1,:),coords(2,:),rows,cols);
        seg = seg + j*double(mask);
    end
    seg = seg(1:mult:end,1:mult:end)+1;
    outfile_name = fullfile(user_annot_dir,[lower(fname(1:end-4)) '.mat']);
    fprintf('Done with image %s in %.2f seconds\n',fname(1:end-4),toc(T));
    parsave(outfile_name,{seg});
end

%{

luong
%}

display('Done');