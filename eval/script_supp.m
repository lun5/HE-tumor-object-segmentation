% this script generate tex for the supplementary entries
SUPP_DIR = '/Users/lun5/Research/data/evaluation_results/supp_materials/';
dir_list = {'best_bdry_300_May30_well_defined','best_bdry_300_May30_invasive',...
    'best_bdry_May31_overlap_well_defined','best_bdry_May31_overlap_invasive'};

fid = fopen('graphic_export.txt','w');
for i = 1: length(dir_list)
   dir_name = fullfile(SUPP_DIR,dir_list{i});
   imlist = dir(fullfile(dir_name,'*.png'));
   imlist = {imlist.name}';
   for j = 1:length(imlist)
       if mod(j,4) == 1
           fprintf(fid,'%s\n%s\n%s\n%s\n','\begin{figure*}','\begin{center}',...
               '\begin{tabular}{|c|}','\hline');
       end
       imname = fullfile(dir_name,imlist{j});
       fprintf(fid,'%s%s%s\n','\includegraphics[width=0.95\textwidth]{',...
           imname,'}\\');
       fprintf(fid,'%s\n','\hline');
       if mod(j,4) == 0
           fprintf(fid,'%s\n%s\n%s\n\n\n','\end{tabular}','\end{center}','\end{figure*}');
       end
       
       if mod(j,20) == 0; fprintf(fid,'%s\n','\clearpage');end
   end    
   fprintf(fid,'\n\n\nDONE\n');
end

fclose(fid);