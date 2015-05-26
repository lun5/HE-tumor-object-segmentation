%% function processInputClassifier(inputFile, which_features)

function processInputClassifier(inputFile, outputDir, which_features)
data = load(inputFile,'sample'); % Read the sample file
sample = data.sample; clear data;
numImages = length(sample);
feat_sample = cell(numImages,1);
opts = setEnvironment_affinity; scale = 1;

fprintf('\nProcessing feature type ''%s'':\n',which_features);
fprintf('Calculate features, differences, and internal statistics...\n');
if strcmp(which_features,'color')
    numFeats = 3;
else
    numFeats = 1 ;
end
numImages = 1; % FOR TESTING
external_stats_sample = cell(numImages,1);
tic;
for ind = 1: numImages
    RGBsamples = sample{ind};
    numSamples = size(RGBsamples,1);
    RGB_A = reshape(RGBsamples(:,3:5),[1 numSamples 3]);
    RGB_B = reshape(RGBsamples(:,6:8),[1 numSamples 3]);
    [f_maps] = getFeatures(double(RGB_A)./255,scale,{which_features},opts);
    feat_A = reshape(f_maps{1},[numSamples numFeats]);
    [f_maps] = getFeatures(double(RGB_B)./255,scale,{which_features},opts);
    feat_B = reshape(f_maps{1},[numSamples numFeats]);
    % calculate distance
    if ~ strcmp(which_features,'hue opp')
        d = sqrt(sum((feat_A-feat_B).^2,2));
    else
        d = circ_dist(feat_A,feat_B);
    end
    F = [feat_A feat_B];
    F_unary = [F(:,1:numFeats); F(:,numFeats+1:2*numFeats)]; 
    A_idx = 1:numSamples; B_idx = (numSamples+1):(2*numSamples);
    % calculate logP(A,B) internal
    if ~ strcmp(which_features,'hue opp')
        [p] = learnP_A_Blocal(F, opts);
        [pmi,pJoint,~] = evalPMI(p,F,F_unary,A_idx,B_idx,opts);
    else
        % calculate pmi internal for hue opp
        [mixture_params] = learnMixtureModel(F,opts);
        [pmi,pJoint,~] = evalPMI_theta(F,mixture_params,opts);
    end
    feat_sample{ind} = [sample{ind}(:,end) feat_A feat_B d log(pJoint) pmi];
    ind_sample_external = randperm(numSamples,200);
    external_stats_sample{ind} = [feat_A(ind_sample_external,:) feat_B(ind_sample_external,:)];
    disp(ind)
end
t = toc; fprintf('done: %1.2f sec\n', t);
clear sample;
external_stats_sample = cat(1,external_stats_sample{1:end});
fprintf('Estimate external statistics...\n');tic;
%% learn external statistics
if ~ strcmp(which_features,'hue opp')
    [p_external] = learnP_A_Blocal(external_stats_sample, opts);
    mixture_params_external = [];
else
    % calculate pmi internal for hue opp
    [mixture_params_external] = learnMixtureModel(external_stats_sample,opts);
    p_external = [];
end
clear external_stats_sample
t = toc; fprintf('done: %1.2f sec\n', t);

fprintf('Calculate external statistics for each image ...\n');tic;
for ind = 1:numImages
    F = feat_sample{ind}(:,2:(2*numFeats+1));
    numSamples = size(F,1);
    F_unary = [F(:,1); F(:,2)]; A_idx = 1:numSamples; 
    B_idx = (numSamples+1):(2*numSamples);
    if ~ strcmp(which_features,'hue opp')
        [pmi_external,pJoint_external,~] = evalPMI(p_external,F,F_unary,A_idx,B_idx,opts);
    else
        % calculate pmi internal for hue opp
        [pmi_external,pJoint_external,~] = evalPMI_theta(F,mixture_params_external,opts);
    end
    feat_sample{ind} = cat(2,feat_sample{ind},log(pJoint_external),pmi_external);
    disp(ind)
end
t = toc; fprintf('done: %1.2f sec\n', t);

splitStr = regexp(which_features,' ','split');
fname = fullfile(outputDir,[splitStr{1} '_cue.mat']);
save(fname,'feat_sample');

end

