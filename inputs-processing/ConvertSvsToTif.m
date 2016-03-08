


noOfImages=0;
imageFolder = '/home/tosun/ImageData/LUNG_fromFSU/8controlFFPE/svs/';



dirName= imageFolder;
dirData = dir( fullfile(dirName,'*.svs') );      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
        fileList,'UniformOutput',false);
end
subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
%#   that are not '.' or '..'
for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(pwd,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; getAllFiles(nextDir)];  %# Recursively call getAllFiles
end
noOfImages =noOfImages +length(fileList);
for j=1:length(fileList)
    system(['vips extract_band ' fileList{j} ' ' fileList{j}(1:end-4) '.tif[pyramid,tile,compression=jpeg,tile-width=256,tile-height=256] 0 --n 3']);
    done_image = fileList{j}
end


