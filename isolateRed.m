function mask = isolateRed(tmpim)

load(fullfile(pwd,'normalization','ClusterCenters_Breast.mat'));
tmpim = double(tmpim);
%tmpim= imread(fileName);
mrows  = size(tmpim,1);
mcols= size(tmpim,2);

red = tmpim(:,:,1); green = tmpim(:,:,2); blue = tmpim(:,:,3);
vecred = repmat(red(:),[1,size(kMeansClusterCenters,1)]);
vecgreen = repmat(green(:),[1,size(kMeansClusterCenters,1)]);
vecblue = repmat(blue(:),[1,size(kMeansClusterCenters,1)]);
distred = (vecred - repmat(kMeansClusterCenters(:,1)',[numel(red),1])).^2;
distgreen = (vecgreen - repmat(kMeansClusterCenters(:,2)',[numel(red),1])).^2;
distblue = (vecblue - repmat(kMeansClusterCenters(:,3)',[numel(red),1])).^2;
distance = sqrt(double(distred+distgreen+distblue));
[~,label_vector] = min(distance,[],2);
pixel_labels = reshape(label_vector,size(red));
% cluster label 11 is red
mask = (pixel_labels == 11);
end