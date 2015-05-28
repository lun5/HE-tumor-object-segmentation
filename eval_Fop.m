%fname = '95f7k8LoesYevi.mat'; %
%fname = 'aNaggwovpxANWq0.mat';
gt_dir = '/Users/lun5/Research/data/gt_dir';
seg_dir_all = '/Users/lun5/Research/data/segmentResults';
%seg_dir = fullfile(seg_dir_all,'rho200sig25');
seg_dir = fullfile(seg_dir_all,'rho300sig25');
fnames = dir(fullfile(gt_dir,'*.mat'));
fnames =  {fnames.name}';
thresh = 0.01:0.01:1;
Fop_measure_thr = cell(length(fnames),1);
Fop_measure_stat = zeros(length(thresh),4);
Fop_measure_best = zeros(1,4);
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
        %figure; imagesc(seg); axis equal; axis off;
        measure = eval_segm( seg, ground_truth, 'fop' );
        im_fop_th(j,:) = cat(2, thresh(j), measure(1:3));
        %fprintf('thresh %1.3f Fop %1.2f %1.2f %1.2f\n',thresh(j),measure(1:3));
    end
    Fop_measure_thr{i} = im_fop_th;
    Fop_measure_stat = Fop_measure_stat + im_fop_th./numImages;
    [~,ind_best] = max(im_fop_th(:,2));
    Fop_measure_best = Fop_measure_best + im_fop_th(ind_best,:)./numImages;
    i
end

[Fmax, ind] = max(Fop_measure_stat(:,2));
Fop_measure_stat(ind,:)
figure; plot(Fop_measure_stat(:,4),Fop_measure_stat(:,3));hold on;
plot(Fop_measure_stat(ind,4), Fop_measure_stat(ind,3),'ro','MarkerSize',15,'MarkerFaceColor','r');hold off;

Fop_measure_best(2) = 2*(Fop_measure_best(3)*Fop_measure_best(4))/(Fop_measure_best(3) + Fop_measure_best(4));
Fop_measure_best

Fop_measure_all = cat(1,Fop_measure_thr{1:end});
[F_max_ois, ind_ois] = max(Fop_measure_all(:,2));
Fop_measure_all(ind_ois,:)

Fop_best = Fop_measure_thr{32};
figure; plot(Fop_best(:,4),Fop_best(:,3));hold on;
plot(Fop_best(mod(ind_ois,100),4), Fop_best(mod(ind_ois,100),3),'ro','MarkerSize',15,'MarkerFaceColor','r');hold off;
