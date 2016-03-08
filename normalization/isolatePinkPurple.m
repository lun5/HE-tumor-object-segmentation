function mask = isolatePinkPurple(tmpim)

load(fullfile(pwd,'normalization','ClusterCenters_Breast.mat'));

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

for i = 1: mrows
    for j = 1: mcols
        if (pixel_labels(i, j) == 10) || (pixel_labels(i, j) == 11) 
            pixel_labels(i, j)= 0;
        else
            pixel_labels(i, j)= 1;
        end
    end
end

mask= pixel_labels;

end