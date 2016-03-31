%% convert non expert annotations to polygon for further processing
IMG_DIR = 'Z:\TilesForLabeling_tiff_renamed';

annot_dir = 'Z:\HEproject\data\GroundTruth\non_expert';
user_name = 'Maurice';
%user_name = 'Om';

user_annot_dir = fullfile(annot_dir,user_name);

fname = fullfile(user_annot_dir,'annotations.txt');
text = fileread(fname);text = strrep(text,'"','');
coords_files = strsplit(text,'polygons'); %Maurice
%coords_files = strsplit(text,'filename:'); %Om
coords_files = coords_files(2:end);
num_files = length(coords_files);
rows = 2048; cols = 2048;
mult = 4;
im_list = {'1bhjqxecct','3tdnf6bhanan','4nkj5wqcqj',...
    '9fukmjujwi','13nedzdzfh','apaal7fc2payi','aqizfuqbbxyu','bylklqnsvf4d',...
    'drfmkoerzy','finkidqlnznihk','hrdqmlu2ig','ovvp7sjwtrrax3t','sfzywlg4291ch9',...
    'shyy81fuot9hu','sp2byg1b33ghl','t42x5lgm4py','taljyo23jlxd','tj6ac55wqwfi',...
    'uhin9nl4ju7bls','updte6rxh7afsf','uraxeh1spli7ky9', 'uyuukdfzqq','vdqmu8xq2xc',...
    'vm3qo9caekfodi','vwxsdt6g3n2f','vyggwhmvuj1','w8kwtop6hyp','w9fpyfgxtibvdk',...
    'wm9g9oyreu','xhcy7eroyqz','ycivjoxn14stvq','ylf5dghxyp2sxp','ynqnt1mljk4e',...
    'yxr1wm9fonpvmf','yxt9szgtneh8','z8wtvu0v8g','zlelrgflyyft4c','zxvcwqjoeyd'};
for i = 1:num_files
    T = tic;
    coordinates = regexp(coords_files{i},'coords:','split');
    fname = regexp(coordinates{end},'filename:','split'); %Maurice
    coordinates{end} = fname{1}; 
    fname = strtrim(regexprep(fname{end},'[},{]',''));
    fname = lower(fname(1:end-4));
    %coordinates = regexp(coords_files{i},'coords:','split');
    %fname = strsplit(coordinates{1},','); % Om
    %fname = strtrim(fname{1}); fname = lower(fname(1:end-5));
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
    if ismember(fname,im_list)
        seg = seg(1:512,1:512)+1;
    else
        seg = seg(1:mult:end,1:mult:end)+1;
    end
    outfile_name = fullfile(user_annot_dir,[fname '.mat']);
    fprintf('Done with image %s in %.2f seconds\n',fname,toc(T));
    parsave(outfile_name,{seg});
end

%{

luong
%}

display('Done');