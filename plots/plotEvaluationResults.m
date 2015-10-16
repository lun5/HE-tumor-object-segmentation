%% Plot the evaluation for luminance, col+var, hue 
% luminance_dir = '/Users/lun5/Research/data/eval_luminance';
% col_var_dir = '/Users/lun5/Research/data/eval_col_var_512_sig3_exp2';
% hue_dir = '/Users/lun5/Research/data/eval_hue_512_512_sig3_exp2_newAffinity_Merge/finalResults';
result_dir = 'Z:\HEproject\evaluation_results';
evaluation_dir = {fullfile(result_dir,'bsr','ev_txt_reannotated_Oct12'),...
    fullfile(result_dir,'Isola_lowres_accurate','ev_txt_reannotated_Oct12'),...
    fullfile(result_dir,'Isola_speedy','ev_txt_reannotated_Oct12'),...
    fullfile(result_dir,'eval_PMI_scale_offset_all3_newsetup','ev_txt'),...
    fullfile(result_dir,'JSEG','one_scale','ev_txt_reannotated_Oct14'),...
    fullfile(result_dir,'JSEG','multi_scale','ev_txt_reannotated_Oct14'),...
    fullfile(result_dir,'ncuts_color','ev_txt_reannotated_Oct14'),...
    fullfile(result_dir,'GraphRLM','ev_txt_reannotated_Oct14'),...
    fullfile(result_dir,'EGB','ev_txt_reannotated_Oct14')};

method_names = {'gPb-OWT-UCM','Isola lowres accurate','Isola speedy','opp color','JSEG one scale',...
   'JSEG multiscale','ncuts multiscale','GraphRLM','EGB'};
% method_names = {'gPb-OWT-UCM','Isola lowres accurate','Isola speedy','opp color'};
num_dir = length(evaluation_dir);
bdry_prvals = cell(num_dir,1);
bdry_evalRes = cell(num_dir,1);
Fop_prvals = cell(num_dir,1);
Fop_evalRes = cell(num_dir,1);
cols = rand(num_dir,3);
for i = 1: num_dir 
    bdry_prvals{i} = dlmread(fullfile(evaluation_dir{i},'eval_bdry_thr.txt'));% thresh,r,p,f
    [~,sort_ind] = sort(bdry_prvals{i}(:,2),'ascend');
    bdry_prvals{i} = bdry_prvals{i}(sort_ind,:);
    bdry_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_bdry.txt'));
    
    Fop_prvals{i} = dlmread(fullfile(evaluation_dir{i},'eval_Fop_thr.txt')); % thresh,r,p,f
    [~,sort_ind] = sort(Fop_prvals{i}(:,4),'ascend');
    Fop_prvals{i} = Fop_prvals{i}(sort_ind,:);
    Fop_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_Fop.txt'));
end

%% plot F_boundary
h = [];open('isoF_new.fig');hold on
for i = 1:num_dir
    h(i) = plot(bdry_prvals{i}(1:end,2),bdry_prvals{i}(1:end,3),'Color',cols(i,:),'LineWidth',3);hold on
    plot(bdry_evalRes{i}(2),bdry_evalRes{i}(3),'o','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',8);
    
end
legend(h, method_names);
hold off;
ax = gca;ax.XTick = 0:0.2:1; ax.YTick = 0.2:0.2:1;
set(gca,'FontSize',16);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')


%% plot F_op
h =[]; open('isoF_new.fig');hold on
for i = 1:num_dir
    %$$witch this after the new F_op
    h(i) = plot(Fop_prvals{i}(1:end,4),Fop_prvals{i}(1:end,3),'Color',cols(i,:),'LineWidth',3);hold on
    plot(Fop_evalRes{i}(2),Fop_evalRes{i}(3),'o','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(1,:),'MarkerSize',8);
end
legend(h, method_names);
hold off;
ax = gca;ax.XTick = 0:0.2:1; ax.YTick = 0.2:0.2:1;
set(gca,'FontSize',16);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')

%% plot F_boundary
luminance_prvals = dlmread(fullfile(luminance_dir,'eval_bdry_thr.txt')); % thresh,r,p,f
col_var_prvals = dlmread(fullfile(col_var_dir,'eval_bdry_thr.txt')); % thresh,r,p,f
hue_prvals = dlmread(fullfile(hue_dir,'eval_bdry_thr.txt')); % thresh,r,p,f

luminance_evalRes = dlmread(fullfile(luminance_dir,'eval_bdry.txt'));
col_var_evalRes = dlmread(fullfile(col_var_dir,'eval_bdry.txt'));
hue_evalRes = dlmread(fullfile(hue_dir,'eval_bdry.txt'));

%col = {'r','g','b'};
col = [1 0 0; 0 1 1; 0 0 1];
open('isoF_new.fig');hold on
h1 = plot(luminance_prvals(1:end,2),luminance_prvals(1:end,3),'Color',col(1,:),'LineWidth',3);hold on
h2 = plot(col_var_prvals(1:end,2),col_var_prvals(1:end,3),'Color',col(2,:),'LineWidth',3);
h3 = plot(hue_prvals(1:end,2),hue_prvals(1:end,3),'Color',col(3,:),'LineWidth',3);
%legend([h1,h2,h3],'luminance','col+var','opponent hue');
plot(luminance_evalRes(2),luminance_evalRes(3),'o','MarkerFaceColor',col(1,:),'MarkerEdgeColor',col(1,:),'MarkerSize',8);
plot(col_var_evalRes(2),col_var_evalRes(3),'o','MarkerFaceColor',col(2,:),'MarkerEdgeColor',col(2,:),'MarkerSize',8);
plot(hue_evalRes(2),hue_evalRes(3),'o','MarkerFaceColor',col(3,:),'MarkerEdgeColor',col(3,:),'MarkerSize',8);
hold off
ax = gca;ax.XTick = 0:0.2:1; ax.YTick = 0.2:0.2:1;
set(gca,'FontSize',30);

results_dir = '/Users/lun5/Dropbox/NIPS2015/images';
set(gcf,'PaperPositionMode','auto')
set(gcf,'color','white');
print(fullfile(results_dir,'F_b'),'-dtiff','-r300');

F_boundaries = [luminance_evalRes(4), luminance_evalRes(7);...
    col_var_evalRes(4), col_var_evalRes(7);...
    hue_evalRes(4), hue_evalRes(7)];

%% plot F_object_part
luminance_prvals = dlmread(fullfile(luminance_dir,'eval_Fop_thr.txt')); % thresh,r,p,f
col_var_prvals = dlmread(fullfile(col_var_dir,'eval_Fop_thr.txt')); % thresh,r,p,f
hue_prvals = dlmread(fullfile(hue_dir,'eval_Fop_thr.txt')); % thresh,r,p,f

luminance_evalRes = dlmread(fullfile(luminance_dir,'eval_Fop.txt'));
col_var_evalRes = dlmread(fullfile(col_var_dir,'eval_Fop.txt'));
hue_evalRes = dlmread(fullfile(hue_dir,'eval_Fop.txt'));

open('isoF_new.fig');hold on
h1 = plot(luminance_prvals(1:end,4),luminance_prvals(1:end,3),'Color',col(1,:),'LineWidth',3);
h2 = plot(col_var_prvals(1:end,4),col_var_prvals(1:end,3),'Color',col(2,:),'LineWidth',3);
h3 = plot(hue_prvals(1:end,4),hue_prvals(1:end,3),'Color',col(3,:),'LineWidth',3);
legend([h1,h2,h3],'luminance','col+var','H&E-Hue','Location','southeast');legend boxoff
plot(luminance_evalRes(2),luminance_evalRes(3),'o','MarkerFaceColor',col(1,:),'MarkerEdgeColor',col(1,:),'MarkerSize',8);
plot(col_var_evalRes(2),col_var_evalRes(3),'o','MarkerFaceColor',col(2,:),'MarkerEdgeColor',col(2,:),'MarkerSize',8);
plot(hue_evalRes(2),hue_evalRes(3),'o','MarkerFaceColor',col(3,:),'MarkerEdgeColor',col(3,:),'MarkerSize',8);
hold off
ax = gca;ax.XTick = 0:0.2:1; ax.YTick = 0.2:0.2:1;
set(gca,'FontSize',30);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')
print(fullfile(results_dir,'F_op'),'-dtiff','-r300');

F_op = [luminance_evalRes(4), luminance_evalRes(7);...
    col_var_evalRes(4), col_var_evalRes(7);...
    hue_evalRes(4), hue_evalRes(7)];

%% Covering
luminance_evalRes = dlmread(fullfile(luminance_dir,'eval_cover.txt'));
col_var_evalRes = dlmread(fullfile(col_var_dir,'eval_cover.txt'));
hue_evalRes = dlmread(fullfile(hue_dir,'eval_cover.txt'));

GT_covering = [luminance_evalRes(2), luminance_evalRes(3);...
    col_var_evalRes(2), col_var_evalRes(3);...
    hue_evalRes(2), hue_evalRes(3)];

%% Covering
luminance_evalRes = dlmread(fullfile(luminance_dir,'eval_RI_VOI.txt'));
col_var_evalRes = dlmread(fullfile(col_var_dir,'eval_RI_VOI.txt'));
hue_evalRes = dlmread(fullfile(hue_dir,'eval_RI_VOI.txt'));

PRI = [luminance_evalRes(2), luminance_evalRes(3);...
    col_var_evalRes(2), col_var_evalRes(3);...
    hue_evalRes(2), hue_evalRes(3)];

VOI = [luminance_evalRes(5), luminance_evalRes(6);...
    col_var_evalRes(5), col_var_evalRes(6);...
    hue_evalRes(5), hue_evalRes(6)];

fprintf('F_b(ODS)\tFb(OIS)\tF_op(ODS)\tF_op(OIS)\tCovering(ODS)\tCovering(OIS)\tPRI(ODS)\t PRI(OIS)\t VOI(ODS)\t VOI(OIS)\n')
[F_boundaries F_op GT_covering PRI VOI]