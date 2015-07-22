%% revert the de-identifying 
% add 5 characters in front of the current name to identify the whole slide
% that the tiles comes from
% 06/30/2015

file_mapping = readtable('Z:\mappingNames.txt');
old_names = file_mapping.Old_Names;
new_names = file_mapping.New_Names;
split_old_names = cellfun(@(x) regexp(x,'_','split'), old_names,'UniformOutput',false);
split_old_names = cat(1,split_old_names{:});
wsi_names = unique(split_old_names(:,1));
%% build a mapping for whole slide images
symbols = ['a':'z'  '0':'9'];
MAX_ST_LENGTH = 5;MIN_ST_LENGTH = 3;
num_wsi = length(wsi_names); 

renamed_wsi = wsi_names;
st = '';
for i = 1:num_wsi
    stLength = randi([MIN_ST_LENGTH,MAX_ST_LENGTH]);
    nums = randi(numel(symbols),[1 stLength]);
    st = symbols(nums);
    while ismember(st,renamed_wsi)
        nums = randi(numel(symbols),[1 stLength]);
        st = symbols(nums);
    end    
    renamed_wsi{i} = st;
end
mapping_wsi = [wsi_names,renamed_wsi];
T = cell2table(mapping_wsi,'VariableNames',{'WSI_Names','Renamed_WSI'});
writetable(T,fullfile('Z:\','mappingWSINames.txt'));

%% remapping the tiles to wsi
num_images = length(new_names);
new_tilenames = new_names;
for i = 1:num_images
   wsi_name = split_old_names(i,1);
   [~,indx_wsi] = ismember(wsi_name,wsi_names); 
   wsi_renamed = renamed_wsi{indx_wsi};
   new_tilenames{i} = [wsi_renamed '_' new_names{i}];
end

mapping_newtiles = [new_names, new_tilenames];

T = cell2table(mapping_newtiles,'VariableNames',{'tile_Names','tile_Names_with_wsi'});
writetable(T,fullfile('Z:\','mappingTilesNames.txt'));
