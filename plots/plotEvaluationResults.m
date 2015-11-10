%% Plot the evaluation for luminance, col+var, hue 
% luminance_dir = '/Users/lun5/Research/data/eval_luminance';
% col_var_dir = '/Users/lun5/Research/data/eval_col_var_512_sig3_exp2';
% hue_dir = '/Users/lun5/Research/data/eval_hue_512_512_sig3_exp2_newAffinity_Merge/finalResults';
result_dir = 'Z:\HEproject\evaluation_results';
%result_dir = '/Users/lun5/Research/data/evaluation_results';
evaluation_dir = {fullfile(result_dir,'eval_PMI_hue_offset','ev_txt_Oct'),...
    fullfile(result_dir,'eval_PJoint_hue_fullscale','ev_txt_Oct'),...
    fullfile(result_dir,'Isola_lowres_accurate','ev_txt_reannotated_Oct12'),...
    fullfile(result_dir,'Isola_speedy','ev_txt_reannotated_Oct12'),...
    fullfile(result_dir,'bsr','ev_txt_reannotated_Oct12'),...  
    fullfile(result_dir,'JSEG','new_params','scale1','ev_txt'),...
    fullfile(result_dir,'JSEG','new_params','scale2','ev_txt'),... 
    fullfile(result_dir,'ncut_multiscale_1_6','ev_txt'),...
    fullfile(result_dir,'EGB','seism_params','ev_txt'),...   
    fullfile(result_dir,'QuadTree','ev_txt'),...   
    fullfile(result_dir,'MeanShift','ev_txt'),...
    fullfile(result_dir,'GraphRLM','new_params','ev_txt')};
%  
% human_dir = {fullfile(result_dir,'eval_non_expert','Maurice','ev_txt'),...
%     fullfile(result_dir,'eval_non_expert','Om','ev_txt')};
% 
% num_dir_human = length(human_dir);
% human_bdry = cell(num_dir_human,1);
% human_Fop = cell(num_dir_human,1);
% human_cover = cell(num_dir_human,1);human_ri_voi = cell(num_dir_human,1);
% human_PRI =  zeros(num_dir,2);human_VOI = zeros(num_dir,2);human_GT_covering = zeros(num_dir,2);
% for i = 1:length(human_dir)
%     human_bdry{i} = dlmread(fullfile(human_dir{i},'eval_bdry.txt'));
%     human_Fop{i} = dlmread(fullfile(human_dir{i},'eval_Fop.txt'));    
%     human_cover{i} = dlmread(fullfile(human_dir{i},'eval_cover.txt'));
%     human_ri_voi{i} = dlmread(fullfile(human_dir{i},'eval_RI_VOI.txt'));
% end
    
% evaluation_dir = {fullfile(result_dir,'eval_PMI_hue_offset','ev_txt_all232'),...
%     fullfile(result_dir,'eval_PJoint_hue_fullscale','ev_txt_all232'),...
%     fullfile(result_dir,'Isola_lowres_accurate','ev_txt_all232'),...
%     fullfile(result_dir,'Isola_speedy','ev_txt_all232'),...  
%     fullfile(result_dir,'bsr','ev_txt_all232'),...
%     fullfile(result_dir,'JSEG','new_params','scale1','ev_txt_all232'),...
%     fullfile(result_dir,'JSEG','new_params','scale2','ev_txt_all232'),...
%     fullfile(result_dir,'ncut_multiscale_1_6','ev_txt_all232'),...
%     fullfile(result_dir,'EGB','seism_params','ev_txt_all232'),...
%     fullfile(result_dir,'QuadTree','ev_txt_all232'),...    
%     fullfile(result_dir,'MeanShift','ev_txt_all232'),...
%     fullfile(result_dir,'GraphRLM','new_params','ev_txt_all232')};

method_names = {'H&E-hue-PMI','H&E-hue-PJoint',...
    'Isola-lowres-acc','Isola-speedy','gPb',...
   'JSEG-ss','JSEG-ms','NCut','EGB','QuadTree','MShift','GraphRLM'};
% method_names = {'H&E-hue-PMI','H&E-hue-PJoint',...
%     'Isola-lowres-acc','Isola-speedy','gPb',...
%    'JSEG-ss','JSEG-ms','NCut','EGB','QuadTree','MShift','GraphRLM',...
%    'non expert 1','non expert 2'};

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
    F_boundaries(i,:) = bdry_evalRes{i}([4 7]);
    
    Fop_prvals{i} = dlmread(fullfile(evaluation_dir{i},'eval_Fop_thr.txt')); % thresh,r,p,f
    %[~,sort_ind] = sort(Fop_prvals{i}(:,2),'ascend');
    [~,sort_ind] = sortrows(Fop_prvals{i},[2 -3]);
    Fop_prvals{i} = Fop_prvals{i}(sort_ind,:);
    Fop_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_Fop.txt'));
    F_op(i,:) = Fop_evalRes{i}([4 7]);
    
    cover_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_cover.txt'));
    ri_voi_evalRes{i} = dlmread(fullfile(evaluation_dir{i},'eval_RI_VOI.txt'));
    GT_covering(i,:) = cover_evalRes{i}(2:3);
    PRI(i,:) = ri_voi_evalRes{i}(2:3);
    VOI(i,:) = ri_voi_evalRes{i}(5:6);
end

%% plot F_boundary
cols = distinguishable_colors(10);
plot_ind = [1:5, 8:10];
h = [];open('isoF_new.fig');hold on
for i = plot_ind
    h(i) = plot(bdry_prvals{i}(1:end,2),bdry_prvals{i}(1:end,3),'Color',cols(i,:),'LineWidth',3);hold on
    plot(bdry_evalRes{i}(2),bdry_evalRes{i}(3),'o','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',8);  hold on;  
end

%symbols = {'*','d','o','^'};
symbols = {'h','<','*','d','>','^'};
col_symb = distinguishable_colors(6);
h(6) = plot(bdry_evalRes{6}(2),bdry_evalRes{6}(3),symbols{1},'MarkerFaceColor',col_symb(1,:),'MarkerEdgeColor',col_symb(1,:),'MarkerSize',10);
h(7) = plot(bdry_evalRes{7}(2),bdry_evalRes{7}(3),symbols{2},'MarkerFaceColor',col_symb(2,:),'MarkerEdgeColor',col_symb(2,:),'MarkerSize',10);
h(11) = plot(bdry_evalRes{11}(2),bdry_evalRes{11}(3),symbols{5},'MarkerFaceColor',col_symb(5,:),'MarkerEdgeColor',col_symb(5,:),'MarkerSize',10);
h(12) = plot(bdry_evalRes{12}(2),bdry_evalRes{12}(3),symbols{6},'MarkerFaceColor',col_symb(6,:),'MarkerEdgeColor',col_symb(6,:),'MarkerSize',10);

%h(13) = plot(human_bdry{1}(2),human_bdry{1}(3),'d','MarkerEdgeColor','r','MarkerSize',15,'LineWidth',3);
%h(14) = plot(human_bdry{2}(2),human_bdry{2}(3),'d','MarkerEdgeColor','k','MarkerSize',15,'LineWidth',3);

legend(h, method_names,'Location','northeastoutside');
hold off;
ax = gca;ax.XTick = 0:0.2:1; ax.YTick = 0.2:0.2:1;
set(gca,'FontSize',20);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto');
savedir = 'C:\Users\luong_nguyen\Box Sync\CVPR2016\Segmentation\Paper\latex\images\fig5';
outname = 'bdry_prcurves_79.png';
%outname = 'bdry_prcurves_232_withhuman.png';
print(fullfile(savedir,outname),'-dpng','-r300');
%% plot F_op
h =[]; open('isoF_new.fig');hold on
plot_ind = [1:7,9];
for i = plot_ind
    h(i) = plot(Fop_prvals{i}(1:end,2),Fop_prvals{i}(1:end,3),'Color',cols(i,:),'LineWidth',3);hold on
    plot(Fop_evalRes{i}(2),Fop_evalRes{i}(3),'o','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',8);  hold on;  
end

h(8) = plot(Fop_evalRes{8}(2),Fop_evalRes{8}(3),symbols{3},'MarkerFaceColor',col_symb(3,:),'MarkerEdgeColor',col_symb(3,:),'MarkerSize',10);
h(10) = plot(Fop_evalRes{10}(2),Fop_evalRes{10}(3),symbols{4},'MarkerFaceColor',col_symb(4,:),'MarkerEdgeColor',col_symb(4,:),'MarkerSize',10);
h(11) = plot(Fop_evalRes{11}(2),Fop_evalRes{11}(3),symbols{5},'MarkerFaceColor',col_symb(5,:),'MarkerEdgeColor',col_symb(5,:),'MarkerSize',10);
h(12) = plot(Fop_evalRes{12}(2),Fop_evalRes{12}(3),symbols{6},'MarkerFaceColor',col_symb(6,:),'MarkerEdgeColor',col_symb(6,:),'MarkerSize',10);

%h(13) = plot(human_Fop{1}(2),human_Fop{1}(3),'d','MarkerEdgeColor','r','MarkerSize',15,'LineWidth',3);
%h(14) = plot(human_Fop{2}(2),human_Fop{2}(3),'d','MarkerEdgeColor','k','MarkerSize',15,'LineWidth',3);

legend(h, method_names,'Location','northeastoutside');
hold off;
ax = gca;ax.XTick = 0:0.2:1; ax.YTick = 0.2:0.2:1;
set(gca,'FontSize',20);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')
outname = 'Fop_prcurves_79.png';
%outname = 'Fop_prcurves_232_withhuman.png';
print(fullfile(savedir,outname),'-dpng','-r300');

%% Summary

%fprintf('&F_b(ODS)&Fb(OIS)& F_op(ODS)& F_op(OIS)& Covering(ODS)&Covering(OIS) &PRI(ODS) &PRI(OIS)& VOI(ODS)& VOI(OIS)\\\\ \n')
summary = [F_boundaries F_op GT_covering PRI VOI];
for i = 1:size(summary,1)
    fprintf('%s & %.2f &%.2f &%.2f &%.2f &%.2f &%.2f& %.2f& %.2f& %.2f& %.2f \\\\ \n',method_names{i}, summary(i,:));
end
disp('Done');