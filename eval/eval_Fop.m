%gtDir = '/Users/lun5/Research/data/groundTruth_512_512';
%inDir = '/Users/lun5/Research/data/segmentResults/eval_hue_512_512';
%

function [ Fop_ods, P_ods, R_ods, bestT, Fop_ois, P_ois, R_ois] = eval_Fop(imgDir, gtDir, inDir, outDir,nthresh)

%fnames = dir(fullfile(gtDir,'*.mat'));
fnames = dir(fullfile(imgDir,'*.tif'));
fnames =  {fnames.name}';
fnames = lower(fnames);
%thresh = 0.01:0.01:1;
numImages = length(fnames);
if nargin<5, nthresh = 99; end

Fop_measure_stat = zeros(nthresh,4);
Fop_measure_img = zeros(numImages,5);

% Fb_measure_stat = zeros(nthresh,4);
% Fb_measure_img = zeros(numImages,5);
% 
% SC_measure_stat = zeros(nthresh,2);
% SC_measure_img = zeros(numImages,3);
% 
% PRI_measure_stat = zeros(nthresh,2);
% PRI_measure_img = zeros(numImages,3);
% 
% VOI_measure_stat = zeros(nthresh,2);
% VOI_measure_img = zeros(numImages,3);

for i = 1:numImages
    T = tic;

    tmp = load(fullfile(inDir,[fnames{i}(1:end-4) '.mat']));
    if iscell(tmp.data)
        segs = tmp.data;
    else
        ucm2 = double(tmp.data); %clear tmp;
    end
    tmp = load(fullfile(gtDir,[lower(fnames{i}(1:end-4)),'.mat']));
    ground_truth = uint8(tmp.groundTruth{1,1}.Segmentation); %clear tmp;
    %figure; imagesc(ground_truth); axis off; axis equal; set(gca,'Position',[0 0 1 1]);
    %splitStr = regexp(fnames{i},'\.','split');
    %filename = fullfile(output_dir,[splitStr{1} '_groundtruth.png']);
    %print(gcf, '-dpng', filename;
    
    if exist('ucm2', 'var'),
        pb = double(ucm2(3:2:end, 3:2:end));
    else 
        pb = [];
    end
    if ~exist('segs', 'var')
        thresh = linspace(1/(nthresh+1),1-1/(nthresh+1),nthresh)';
    else
        if nthresh ~= numel(segs)
            warning('Setting nthresh to number of segmentations');
            nthresh = numel(segs);            
        end
        thresh = 1:nthresh; thresh=thresh';
    end
    
    im_fop_th = zeros(nthresh,4);% objects and parts
%     im_fb_th = zeros(nthresh,4);% boundaries
%     im_sc_th = zeros(nthresh, 2); % region covering
%     im_pri_th = zeros(nthresh, 2); % probability rand index
%     im_voi_th = zeros(nthresh, 2); % variation of information
        
    for j = 1:length(thresh)
        if exist('segs','var')
            seg = uint8(segs{j});
        else
            seg = bwlabel(pb <=  thresh(j)); % should I have flip here? Yes
        end
        % object and parts
        measure = eval_segm( seg, ground_truth, 'fop' );
        im_fop_th(j,:) = cat(2, thresh(j), measure(1:3));
%         % boundaries
%         measure = eval_segm( seg, ground_truth, 'fb' );
%         im_fb_th(j,:) = cat(2, thresh(j), measure(1:3));
%         % regions
%         measure = eval_segm( seg, ground_truth, 'sc' );
%         im_sc_th(j,:) = cat(2, thresh(j), measure);
%         
%         measure = eval_segm( seg, ground_truth, 'pri' );
%         im_pri_th(j,:) = cat(2, thresh(j), measure);
%         
%         measure = eval_segm( seg, ground_truth, 'voi' );
%         im_voi_th(j,:) = cat(2, thresh(j), measure);
        
        %fprintf('thresh %1.3f Fop %1.2f %1.2f %1.2f\n',thresh(j),measure(1:3));
    end
    % object and parts
    Fop_measure_stat = Fop_measure_stat + im_fop_th./numImages;
    [bestT,bestR,bestP,bestF] = maxF(thresh,im_fop_th(:,4),im_fop_th(:,3));
    Fop_measure_img(i,:) = cat(2,i,bestT, bestR, bestP, bestF);
%     % boundaries
%     Fb_measure_stat = Fb_measure_stat + im_fb_th./numImages;
%     [bestT,bestR,bestP,bestF] = maxF(thresh,im_fb_th(:,4),im_fb_th(:,3));
%     Fb_measure_img(i,:) = cat(2,i,bestT, bestF, bestP, bestR);
%     % SC
%     SC_measure_stat = SC_measure_stat + im_sc_th./numImages;
%     [~,ind_best] = max(im_sc_th(:,2));
%     SC_measure_img(i,:) = cat(2,i,im_sc_th(ind_best,:));
%     % PRI
%     PRI_measure_stat = PRI_measure_stat + im_pri_th./numImages;
%     [~,ind_best] = max(im_pri_th(:,2));
%     PRI_measure_img(i,:) = cat(2,i,im_pri_th(ind_best,:));
%     % SC
%     VOI_measure_stat = VOI_measure_stat + im_voi_th./numImages;
%     [~,ind_best] = max(im_voi_th(:,2));
%     VOI_measure_img(i,:) = cat(2,i,im_voi_th(ind_best,:));
%    
    %labels2 = bwlabel(ucm <= thresh(ind_best));
    %labels2 = bwlabel(ucm <= 0.2);
    %seg = labels2(2:2:end, 2:2:end);
    %figure; imagesc(seg); axis equal; axis off;set(gca,'Position',[0 0 1 1]);
    %filename = fullfile(output_dir,[splitStr{1} '_bestSeg.png']);
    %print(gcf, '-dpng', filename);close all;
    t = toc(T);
    fprintf('Done with image %s in %.2f seconds\n',fnames{i}(1:end-4),t);
end

%% Object and parts
Fop_measure_stat(:,2) = fmeasure(Fop_measure_stat(:,4),Fop_measure_stat(:,3));
[bestT,R_ods,P_ods,Fop_ods] = maxF(thresh,Fop_measure_stat(:,4),Fop_measure_stat(:,3));

Fop_measure_ois = mean(Fop_measure_img,1);
P_ois = Fop_measure_ois(4); R_ois = Fop_measure_ois(5);
Fop_ois = fmeasure(R_ois,P_ois);

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
dlmwrite(fname, Fop_measure_stat(:,[4 3 2]),'delimiter','\t','precision',3);

% %% boundaries
% Fb_measure_stat(:,2) = fmeasure(Fb_measure_stat(:,4), Fb_measure_stat(:,3));
% [bestT,R_ods,P_ods,Fb_ods] = maxF(thresh,Fb_measure_stat(:,4),Fb_measure_stat(:,3));
% 
% Fb_measure_ois = mean(Fb_measure_img,1);
% P_ois = Fb_measure_ois(4); R_ois = Fb_measure_ois(5);
% Fb_ois = fmeasure(R_ois,P_ois);
% 
% fname = fullfile(outDir,'eval_bdry.txt');
% fid = fopen(fname,'w');
% if fid==-1,
%    error('Could not open file %s for writing.',fname);
% end
% fprintf(fid,'%10g %10g %10g %10g %10g %10g %10g %10g\n',bestT,R_ods,P_ods,Fb_ods,R_ois,P_ois,Fb_ois);
% fclose(fid);
% 
% fname = fullfile(outDir,'eval_bdry_img.txt');
% fid = fopen(fname,'w');
% if fid==-1,
%     error('Could not open file %s for writing.',fname);
% end
% fprintf(fid,'%10d %10g %10g %10g %10g\n',Fb_measure_img');
% fclose(fid);
% 
% fname = fullfile(outDir,'eval_bdry_thr.txt');
% dlmwrite(fname, Fb_measure_stat,'delimiter','\t','precision',3);
% 
% %% regions
% fname = fullfile(outDir,'eval_cover_th.txt');
% fid = fopen(fname,'w');
% if fid==-1,
%     error('Could not open file %s for writing.',fname);
% end
% fprintf(fid,'%10d %10g\n',SC_measure_stat);
% fclose(fid);
% %dlmwrite(fname, SC_measure_stat,'delimiter','\t','precision',3);
% 
% fname = fullfile(outDir,'eval_cover_img.txt');
% fid = fopen(fname,'w');
% if fid==-1,
%     error('Could not open file %s for writing.',fname);
% end
% fprintf(fid,'%10g %10g %10g\n',SC_measure_img');
% fclose(fid);
% %dlmwrite(fname, SC_measure_img,'delimiter','\t','precision',3);
% 
% [bestR,bestT] = max(SC_measure_stat(:,2));
% R_best = mean(SC_measure_img(:,3));
% R_best_total = max(SC_measure_img(:,3));
% 
% fname = fullfile(outDir,'eval_cover.txt');
% fid = fopen(fname,'w');
% if fid==-1,
%     error('Could not open file %s for writing.',fname);
% end
% fprintf(fid,'%10g %10g %10g %10g\n',bestT, bestR, R_best, R_best_total);
% fclose(fid);
% 
% %% PR and VOI
% fname = fullfile(outDir,'eval_RI_VOI_th.txt');
% fid = fopen(fname,'w');
% if fid==-1,
%     error('Could not open file %s for writing.',fname);
% end
% fprintf(fid,'%10g %10g %10g\n',[thresh PRI_measure_stat(:,2) VOI_measure_stat(:,2)]');
% fclose(fid);
% 
% [RI_best, igRI] = max(PRI_measure_img(:,3));
% [VOI_best, igVOI] = max(VOI_measure_img(:,3));
% bgRI = mean(PRI_measure_img(:,3));
% bgVOI = mean(VOI_measure_img(:,3));
% fname = fullfile(outDir,'eval_RI_VOI.txt');
% fid = fopen(fname,'w');
% if fid==-1,
%     error('Could not open file %s for writing.',fname);
% end
% fprintf(fid,'%10g %10g %10g %10g %10g %10g\n',thresh(igRI), bgRI, RI_best, thresh(igVOI),  bgVOI, VOI_best);
% fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute f-measure fromm recall and precision
function [f] = fmeasure(r,p)
f = 2*p.*r./(p+r+((p+r)==0));
end

% interpolate to find best F and coordinates thereof
function [bestT,bestR,bestP,bestF] = maxF(thresh,R,P)
bestT = thresh(1);
bestR = R(1);
bestP = P(1);
bestF = fmeasure(R(1),P(1));
for i = 2:numel(thresh),
  for d = linspace(0,1),
    t = thresh(i)*d + thresh(i-1)*(1-d);
    r = R(i)*d + R(i-1)*(1-d);
    p = P(i)*d + P(i-1)*(1-d);
    f = fmeasure(r,p);
    if f > bestF,
      bestT = t;
      bestR = r;
      bestP = p;
      bestF = f;
    end
  end
end

end
