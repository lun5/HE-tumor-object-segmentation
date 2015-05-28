%fname = '95f7k8LoesYevi.mat'; %
%fname = 'aNaggwovpxANWq0.mat';
gt_dir = fullfile(pwd,'gt_dir');
seg_dir = fullfile(pwd,'segmentResults','rho200sig100');
fnames = dir(fullfile(gt_dir,'*.mat'));
fnames =  {fnames.name}';
thresh = 0.01:0.01:1;
Fop_measure_thr = cell(length(fnames),1);
Fop_measure_stat = zeros(length(thresh),4);
numImages = length(fnames);
parfor i = 1:numImages
    tmp = load(fullfile(seg_dir,fnames{i})); 
    ucm = double(tmp.data); %clear tmp;
    tmp = load(fullfile(gt_dir,fnames{i}));
    ground_truth = tmp.groundTruth{1,1}.Segmentation; %clear tmp;
    im_fop_th = zeros(length(thresh),4);
    for j = 1:length(thresh)
        labels2 = bwlabel(ucm <= thresh(j));
        seg = labels2(2:2:end, 2:2:end);
        figure; imagesc(seg); axis equal; axis off;
        measure = eval_segm( seg, ground_truth, 'fop' );
        im_fop_th(j,:) = cat(2, thresh(j), measure(1:3));
        %fprintf('thresh %1.3f Fop %1.2f %1.2f %1.2f\n',thresh(j),measure(1:3));
    end
    Fop_measure_thr{i} = im_fop_th;
    Fop_measure_stat = Fop_measure_stat + im_fop_th./numImages;
end




