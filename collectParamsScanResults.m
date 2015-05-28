% collect data from scanning parameter
DATA_DIR ='/home/lun5/HEproject/'; % linux
params_scan_dir = fullfile(DATA_DIR,'evaluation_results','eval_hue_scan');
dir_list = dir(fullfile(params_scan_dir,'rho*'));
dir_list = {dir_list.name}';
numDirs = length(dir_list);
params_scan_stat = cell(numDirs,1);

parfor i = 1:numDirs
    dirname = dir_list{i};
    splitStr = regexp(dirname,'rho','split');
    splitStr = regexp(splitStr{2},'sig','split');
    rho = str2double(splitStr{1})/100; 
    sigma = str2double(splitStr{2})/100;
    if exist(fullfile(params_scan_dir,dirname,'eval_bdry.txt'),'file'),
        evalRes = dlmread(fullfile(params_scan_dir,dirname,'eval_bdry.txt'));
        %Fb_ODS = evalRes(4); Fb_OIS = evalRes(7); AP = evalRes(8);
    else
        error(['missing eval_bdry.txt file for ',dirname]);
    end
    if exist(fullfile(params_scan_dir,dirname,'eval_cover.txt'),'file')
        evalRes_cover = dlmread(fullfile(params_scan_dir,dirname,'eval_cover.txt'));
        %SC_ODS = evalRes_cover(2); SC_OIS = evalRes_cover(3);
    else
        error(['missing eval_cover.txt file for ',dirname]);
    end
    if exist(fullfile(params_scan_dir,dirname,'eval_RI_VOI.txt'),'file')
        evalRes_ri = dlmread(fullfile(params_scan_dir,dirname,'eval_RI_VOI.txt'));
        %PRI_ODS = evalRes_ri(2); PRI_OIS = evalRes_ri(3);
        %VOI_ODS = evalRes_ri(5); VOI_OIS = evalRes_ri(6);
    else
        error(['missing eval_RI_VOI.txt for ', dirname]);
    end
    %params_scan_stat{i} = cat(1,rho,sigma,Fb_ODS, Fb_OIS, AP, SC_ODS, SC_OIS, ...
    %    PRI_ODS, PRI_OIS, VOI_ODS, VOI_OIS);
    params_scan_stat{i} = cat(1,rho, sigma, evalRes(4), evalRes(7), evalRes(8),...
        evalRes_cover(2:3)',evalRes_ri(2:3)',evalRes_ri(5:6)');
    fprintf('Done with %s\n',dirname);
end
params_scan_stat = cat(2, params_scan_stat{1:end});
params_scan_stat = sortrows(params_scan_stat')'; % sort according to rho and sigma
[max_Fb_ods, ind_fb_ods] = max(params_scan_stat(3,:));
[max_ap, ind_ap] = max(params_scan_stat(5,:));
[max_sc_ods, ind_sc_ods] = max(params_scan_stat(6,:));
[max_pri_ods,ind_pri_ods] = max(params_scan_stat(8,:));
[max_voi_ods,ind_voi_ods] = min(params_scan_stat(10,:));

Fb_ods = reshape(params_scan_stat(3,:), [8 5]); %row = sigma, col = rho
AP = reshape(params_scan_stat(5,:), [8 5]);
sc_ods = reshape(params_scan_stat(6,:), [8 5]);
pri_ods = reshape(params_scan_stat(8,:), [8 5]);
voi_ods = reshape(params_scan_stat(10,:), [8 5]);
xlabel('\rho'); ylabel('\sigma');

figure; imagesc(AP); title('AP'); colorbar;
ax = gca;ax.XTick = 1:1:5;ax.YTick = 1:1:8;
ax.XTickLabel = {'1','1.25','1.5','2','3'};
ax.YTickLabel =  {'0.25','0.5','1','2','3','5','10','15'};
xlabel('\rho'); ylabel('\sigma');

figure; imagesc(Fb_ods);title('Fb ODS');colorbar;
ax = gca;ax.XTick = 1:1:5;ax.YTick = 1:1:8;
ax.XTickLabel = {'1','1.25','1.5','2','3'};
ax.YTickLabel =  {'0.25','0.5','1','2','3','5','10','15'};
xlabel('\rho'); ylabel('\sigma');

figure; imagesc(sc_ods); title('Seg covering ODS');colorbar;
ax = gca;ax.XTick = 1:1:5;ax.YTick = 1:1:8;
ax.XTickLabel = {'1','1.25','1.5','2','3'};
ax.YTickLabel =  {'0.25','0.5','1','2','3','5','10','15'};
xlabel('\rho'); ylabel('\sigma');

figure; imagesc(pri_ods); title('PRI');colorbar;
ax = gca;ax.XTick = 1:1:5;ax.YTick = 1:1:8;
ax.XTickLabel = {'1','1.25','1.5','2','3'};
ax.YTickLabel =  {'0.25','0.5','1','2','3','5','10','15'};
xlabel('\rho'); ylabel('\sigma');

figure; imagesc(voi_ods); title('VOI');colorbar;
ax = gca;ax.XTick = 1:1:5;ax.YTick = 1:1:8;
ax.XTickLabel = {'1','1.25','1.5','2','3'};
ax.YTickLabel =  {'0.25','0.5','1','2','3','5','10','15'};
xlabel('\rho'); ylabel('\sigma');
