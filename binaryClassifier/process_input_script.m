%% process script
outputDir = 'Z:\HEproject\data\pixelPairs\';
%which_features = 'luminance';
%which_features = 'hue opp';
%which_features = 'color';
feature_set = {'luminance','hue opp','color'};
opts_affinity = setEnvironment_affinity;
opts_affinity.affinity.plot = false;
for i  = 1:length(feature_set)
    which_features = feature_set{i};
    inputFile = fullfile(outputDir,'RGBsample_pairs_labels_sig025.mat');
    processInputClassifier(inputFile, outputDir, which_features);
end

