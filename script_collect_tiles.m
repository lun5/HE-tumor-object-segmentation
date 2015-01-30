% script to collect data tiles
sourcedir = 'Z:\';
svs_fnames = dir(fullfile(sourcedir,'svs','*.svs'));
svs_fnames = {svs_fnames.name}';
num_svs = length(svs_fnames);
% folder storing the tiles of interest
destination_dir = fullfile(sourcedir,'TilesForLabeling');
outfile = fopen('tiles_fnames.txt','w'); % record the tiles' names
if ~exist(destination_dir,'dir')
    mkdir(destination_dir);
    fileattrib(destination_dir,'+w');
end
numtiles = 10;
bubble_svs = {'tp09-16-1.svs','tp09-39-1.svs','tp09-777-1.svs',...
    'tp09-1813-1.svs','tp10-420-1.svs','tp10-420-2.svs'};
for i = 1:num_svs
    svs_fname = svs_fnames{i};
    if sum(ismember(bubble_svs,svs_fname)) > 0
        continue;
    end
    
    tilenames = wsi_get_training( sourcedir, svs_fname, destination_dir, numtiles);
    fprintf(outfile,'%s',tilenames{1,1:end-1});
    fprintf(outfile,'%s\n',tilenames{1,end});
end

fclose(outputfile);

