%% Plot the evaluation for luminance, col+var, hue 
% luminance_dir = '/Users/lun5/Research/data/eval_luminance';
% col_var_dir = '/Users/lun5/Research/data/eval_col_var_512_sig3_exp2';
% hue_dir = '/Users/lun5/Research/data/eval_hue_512_512_sig3_exp2_newAffinity_Merge/finalResults';
%result_dir = 'Z:\HEproject\evaluation_results';
result_dir = '/Users/lun5/Research/data/evaluation_results';
evaluation_dir = {fullfile(result_dir,'opp_col_lowres_PMI_3channels','ev_txt'),...
    fullfile(result_dir,'bsr','ev_txt_reannotated_Oct12'),...
    fullfile(result_dir,'Isola_lowres_accurate','ev_txt_reannotated_Oct12'),...
    fullfile(result_dir,'Isola_speedy','ev_txt_reannotated_Oct12'),...    
    fullfile(result_dir,'JSEG','new_params','scale1','ev_txt'),...
    fullfile(result_dir,'JSEG','new_params','scale2','ev_txt'),... 
    fullfile(result_dir,'ncut_multiscale_1_6','ev_txt'),...
    fullfile(result_dir,'GraphRLM','new_params','ev_txt'),...
    fullfile(result_dir,'EGB','seism_params','ev_txt'),...
    fullfile(result_dir,'QuadTree','ev_txt'),...
    fullfile(result_dir,'MeanShift','ev_txt')};

method_names = {'opp color','gPb-OWT-UCM','Isola lowres accurate','Isola speedy','JSEG one scale',...
   'JSEG multiscale','ncuts multiscale','GraphRLM','EGB','QuadTree','MeanShift'};
%method_names = {'ncuts multiscale','GraphRLM','EGB','QuadTree','MeanShift'};
num_dir = length(evaluation_dir);
% Measures: boundaries, object-parts, region covering, PRI, VOI
bdry_prvals = cell(num_dir,1);
bdry_evalRes = cell(num_dir,1); F_boundaries = zeros(num_dir,2);
Fop_prvals = cell(num_dir,1);
Fop_evalRes = cell(num_dir,1); F_op = zeros(num_dir,2);
cover_evalRes = cell(num_dir,1);ri_voi_evalRes = cell(num_dir,1);
PRI =  zeros(num_dir,2);VOI = zeros(num_dir,2);GT_covering = zeros(num_dir,2);

for i = 1: num_dir 
    bdry_prvals{i} = dlmread(fullfile(evaluation_dir{i},'eval_bdry_thr.txt'));% thresh,r,p,f
    %[~,sort_ind] = sort(bdry_prvals{i}(:,2),'ascend');
    [~,sort_ind] = sortrows(bdry_prvals{i},[2 -3]);
    bdry_prvals{i} = bdry_prvals{i}(sort_ind,:);
    bdry_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_bdry.txt'));
    F_boundaries(i,:) = bdry_evalRes{i}([4 6]);
    
    Fop_prvals{i} = dlmread(fullfile(evaluation_dir{i},'eval_Fop_thr.txt')); % thresh,r,p,f
    %[~,sort_ind] = sort(Fop_prvals{i}(:,2),'ascend');
    [~,sort_ind] = sortrows(Fop_prvals{i},[2 -3]);
    Fop_prvals{i} = Fop_prvals{i}(sort_ind,:);
    Fop_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_Fop.txt'));
    F_op(i,:) = Fop_evalRes{i}([4 6]);
    
    cover_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_cover.txt'));
    ri_voi_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_RI_VOI.txt'));
    GT_covering(i,:) = cover_evalRes{i}(2:3);
    PRI(i,:) = ri_voi_evalRes{i}(2:3);
    VOI(i,:) = ri_voi_evalRes{i}(5:6);
end

%% plot F_boundary
cols = rand(num_dir,3);
plot_ind = 8:num_dir;
h = [];open('isoF_new.fig');hold on
for i = plot_ind
    h(i) = plot(bdry_prvals{i}(1:end,2),bdry_prvals{i}(1:end,3),'Color',cols(i,:),'LineWidth',3);hold on
    plot(bdry_evalRes{i}(2),bdry_evalRes{i}(3),'o','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',8);
    
end
legend(h(plot_ind), method_names(plot_ind));
hold off;
ax = gca;ax.XTick = 0:0.2:1; ax.YTick = 0.2:0.2:1;
set(gca,'FontSize',16);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')


%% plot F_op
h =[]; open('isoF_new.fig');hold on
plot_ind = 5:num_dir;
for i = plot_ind
    %$$witch this after the new F_op
    %h(i) = plot(Fop_prvals{i}(1:end,4),Fop_prvals{i}(1:end,3),'Color',cols(i,:),'LineWidth',3);hold on
    h(i) = plot(Fop_prvals{i}(1:end,2),Fop_prvals{i}(1:end,3),'Color',cols(i,:),'LineWidth',3);hold on
    plot(Fop_evalRes{i}(2),Fop_evalRes{i}(3),'o','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',8);
end
legend(h(plot_ind), method_names(plot_ind));
hold off;
ax = gca;ax.XTick = 0:0.2:1; ax.YTick = 0.2:0.2:1;
set(gca,'FontSize',16);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')

%% Summary

fprintf('F_b(ODS)\tFb(OIS)\tF_op(ODS)\tF_op(OIS)\tCovering(ODS)\tCovering(OIS)\tPRI(ODS)\t PRI(OIS)\t VOI(ODS)\t VOI(OIS)\n')
[F_boundaries F_op GT_covering PRI VOI]