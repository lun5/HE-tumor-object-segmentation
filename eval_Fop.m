%gtDir = '/Users/lun5/Research/data/groundTruth_512_512';
%inDir = '/Users/lun5/Research/data/segmentResults/eval_hue_512_512';
%

function [ Fop_ods, P_ods, R_ods, bestT, Fop_ois, P_ois, R_ois] = eval_Fop(gtDir, inDir, outDir)

fnames = dir(fullfile(gtDir,'*.mat'));
fnames =  {fnames.name}';
thresh = 0.01:0.01:1;
%Fop_measure_thr = cell(length(fnames),1);
Fop_measure_stat = zeros(length(thresh),4);
numImages = length(fnames);
Fop_measure_img = zeros(numImages,5);
parfor i = 1:numImages
    tmp = load(fullfile(inDir,fnames{i})); 
    ucm = double(tmp.data); %clear tmp;
    tmp = load(fullfile(gtDir,fnames{i}));
    ground_truth = tmp.groundTruth{1,1}.Segmentation; %clear tmp;
    %figure; imagesc(ground_truth); axis off; axis equal; set(gca,'Position',[0 0 1 1]);
    %splitStr = regexp(fnames{i},'\.','split');
    %filename = fullfile(output_dir,[splitStr{1} '_groundtruth.png']);
    %print(gcf, '-dpng', filename);
    im_fop_th = zeros(length(thresh),4);
    for j = 1:length(thresh)
        labels2 = bwlabel(ucm <= thresh(j));
        seg = labels2(2:2:end, 2:2:end);
        %figure; imagesc(seg); axis equal; axis off;
        measure = eval_segm( seg, ground_truth, 'fop' );
        im_fop_th(j,:) = cat(2, thresh(j), measure(1:3));
        %fprintf('thresh %1.3f Fop %1.2f %1.2f %1.2f\n',thresh(j),measure(1:3));
    end
    %Fop_measure_thr{i} = im_fop_th;
    Fop_measure_stat = Fop_measure_stat + im_fop_th./numImages;
    [~,ind_best] = max(im_fop_th(:,2));
    %labels2 = bwlabel(ucm <= thresh(ind_best));
    %labels2 = bwlabel(ucm <= 0.2);
    %seg = labels2(2:2:end, 2:2:end);
    %figure; imagesc(seg); axis equal; axis off;set(gca,'Position',[0 0 1 1]);
    %filename = fullfile(output_dir,[splitStr{1} '_bestSeg.png']);
    %print(gcf, '-dpng', filename);close all;
    Fop_measure_img(i,:) = cat(2,i,im_fop_th(ind_best,:));
    i
end

Fop_measure_stat(:,2) = 2*(Fop_measure_stat(:,3).*Fop_measure_stat(:,4))./(Fop_measure_stat(:,3)+Fop_measure_stat(:,4));
[Fop_ods, ind] = max(Fop_measure_stat(:,2));
P_ods = Fop_measure_stat(ind,3); R_ods = Fop_measure_stat(ind,4); bestT = Fop_measure_stat(ind,1);

Fop_measure_ois = mean(Fop_measure_img,1);
Fop_ois = 2*(Fop_measure_ois(4)*Fop_measure_ois(5))/(Fop_measure_ois(4) + Fop_measure_ois(5));
P_ois = Fop_measure_ois(4);R_ois = Fop_measure_ois(5);
fname = fullfile(outDir,'eval_Fop.txt');
fid = fopen(fname,'w');
if fid==-1,
   error('Could not open file %s for writing.',fname);
end
fprintf(fid,'%10g %10g %10g %10g %10g %10g %10g %10g\n',bestT,R_ods,P_ods,Fop_ods,R_ois,P_ois,Fop_ois);
fclose(fid);

fname = fullfile(outDir,'eval_Fop_img.txt');
fid = fopen(fname,'w');
if fid==-1,
    error('Could not open file %s for writing.',fname);
end
fprintf(fid,'%10d %10g %10g %10g %10g\n',Fop_measure_img');
fclose(fid);

fname = fullfile(outDir,'eval_Fop_thr.txt');
dlmwrite(fname, Fop_measure_stat,'delimiter','\t','precision',3);


%open('isoF.fig'); hold on;
%plot(Fop_measure_stat(:,4),Fop_measure_stat(:,3),'LineWidth',3);%hold on;
%plot(Fop_measure_stat(ind,4), Fop_measure_stat(ind,3),'ro','MarkerSize',15,'MarkerFaceColor','r');hold off;
%plot_eval(outDir);

%Fop_measure_all = cat(1,Fop_measure_thr{1:end});
%[F_max_ois, ind_ois] = max(Fop_measure_all(:,2));
%Fop_measure_all(ind_ois,:)

%Fop_best = Fop_measure_thr{32};
%figure; plot(Fop_best(:,4),Fop_best(:,3));hold on;
%plot(Fop_best(mod(ind_ois,100),4), Fop_best(mod(ind_ois,100),3),'ro','MarkerSize',15,'MarkerFaceColor','r');hold off;
