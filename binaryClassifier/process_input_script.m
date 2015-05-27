%% process script
outputDir = 'Z:\HEproject\data\pixelPairs\';
%which_features = 'luminance';
%which_features = 'hue opp';
%which_features = 'color';
feature_set = {'luminance','hue opp','color'};

for i  = 1:length(feature_set)
    which_features = feature_set{i};
    inputFile = fullfile(outputDir,'RGBsample_pairs_labels.mat');
    processInputClassifier(inputFile, outputDir, which_features);
end

