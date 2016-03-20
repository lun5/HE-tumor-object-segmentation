function plot_eval(evalDir,col)
% plot evaluation results.
% Pablo Arbelaez <arbelaez@eecs.berkeley.edu>

if nargin<2, col = 'r'; end

fwrite(2,sprintf('\n%s\n',evalDir));

if exist(fullfile(evalDir,'eval_bdry_thr.txt'),'file'),
    open('isoF_new.fig');
    hold on
    prvals = dlmread(fullfile(evalDir,'eval_bdry_thr.txt')); % thresh,r,p,f
    f=find(prvals(:,2)>=0.01);
    prvals = prvals(f,:);


    evalRes = dlmread(fullfile(evalDir,'eval_bdry.txt'));
    %if size(prvals,1)>1,
        plot(prvals(1:end,2),prvals(1:end,3),col,'LineWidth',3);
    %else
        plot(evalRes(2),evalRes(3),'o','MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',8);
    %end
    hold off
    
    fout = fopen(fullfile(evalDir,'finalResults.txt'),'w')
    fprintf(fout,'Boundary\n');

    fprintf(fout,'ODS: F( %1.2f, %1.2f ) = %1.2f   [th = %1.2f]\n',evalRes(2:4),evalRes(1));
    fprintf(fout,'OIS: F( %1.2f, %1.2f ) = %1.2f\n',evalRes(5:7));
    fprintf(fout,'Area_PR = %1.2f\n\n',evalRes(8));
end

if exist(fullfile(evalDir,'eval_Fop_thr.txt'),'file'),
    open('isoF_new.fig');
    hold on
    prvals = dlmread(fullfile(evalDir,'eval_Fop_thr.txt')); % thresh,r,p,f
    f=find(prvals(:,2)>=0.01);
    prvals = prvals(f,:);
    evalRes = dlmread(fullfile(evalDir,'eval_Fop.txt'));
    %if size(prvals,1)>1,
        plot(prvals(:,2),prvals(:,3),col,'LineWidth',3);
    %else
        plot(evalRes(2),evalRes(3),'o','MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',10);
    %end
    hold off
end

if exist(fullfile(evalDir,'eval_Fop.txt'),'file'),
    evalRes = dlmread(fullfile(evalDir,'eval_Fop.txt'));
    fprintf(fout,'Objects and parts\n');
    fprintf(fout,'ODS: F( %1.2f, %1.2f ) = %1.2f   [th = %1.2f]\n',evalRes(2:4),evalRes(1));

    fprintf(fout,'OIS: F( %1.2f, %1.2f ) = %1.2f\n',evalRes(5:7));
    %fprintf('Area_PR = %1.2f\n\n',evalRes(8));
end

if exist(fullfile(evalDir,'eval_cover.txt'),'file'),
    evalRes = dlmread(fullfile(evalDir,'eval_cover.txt'));
    fprintf(fout,'Region\n');
    fprintf(fout,'GT covering: ODS = %1.2f [th = %1.2f]. OIS = %1.2f. Best = %1.2f\n',evalRes(2),evalRes(1),evalRes(3:4));
    evalRes = dlmread(fullfile(evalDir,'eval_RI_VOI.txt'));
    fprintf(fout,'Rand Index: ODS = %1.2f [th = %1.2f]. OIS = %1.2f.\n',evalRes(2),evalRes(1),evalRes(3));
    fprintf(fout,'Var. Info.: ODS = %1.2f [th = %1.2f]. OIS = %1.2f.\n',evalRes(5),evalRes(4),evalRes(6));

end
