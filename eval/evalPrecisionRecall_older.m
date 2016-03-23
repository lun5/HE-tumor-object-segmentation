
function [precision,recall, penalty] = evalPrecisionRecall_old(groundTruth,result)

% if the segmentation result is just one segment
if max(result(:)) <= 1
    precision = 0;
    recall = 0;
    penalty = 0;
    return;
end
gto = groundTruth{1,1};
bmap = gto.Boundaries;
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
logical_cells = cellfun(cellfind('stroma'),gto.names);
stromaLabel = find(logical_cells==1);
logical_cells = cellfun(cellfind('white'),gto.names);
whiteLabel = find(logical_cells==1);
gt = gto.Segmentation;
noback = gt~=stromaLabel &  gt~=whiteLabel;
logical_cells = cellfun(cellfind('white space'),gto.names);
whiteLabel = find(logical_cells==1);

gt = gto.Segmentation;
gt(bmap) = 0;
if ~isempty(stromaLabel)
    gt(gt==0)=stromaLabel;
end
if ~isempty(whiteLabel) && ~isempty(stromaLabel)
    noback = gt~=stromaLabel &  gt~=whiteLabel;
elseif ~isempty(stromaLabel);
    noback = gt~=stromaLabel;
else 
    noback = ones(size(gt));
end
% WHAT HAPPENS IF THERE IS ONLY STROMA
L= bwlabeln(noback); % NEED TO PUT 0 FOR BOUNDARIES

totalRegion=zeros(max(max(result)),max(max(L)));
for i = 1:size(result,1)
    for j = 1:size(result,2)
        rno= result(i,j);
        gs = L(i,j);
        if gs~=0 && rno~=0
            totalRegion(rno,gs)=totalRegion(rno,gs)+1;
        end
    end
end

totalhit =0;
totalseg=0;
totalorg = 0;
totalhitGT = 0;
totalhitSEG=0;
for i=1:max(max(L))
    thr = sum(sum(L==i))/20;
    indf = find(totalRegion(:,i)>thr);
    if indf>0
        totalhitGT = totalhitGT+1;
    end
    for j=1:length(indf)
%         [hit,in] = max(totalRegion(:,i));
        hit = totalRegion(indf(j),i);
        totalhit = totalhit +hit;
        seg = sum(sum(result==indf(j)));
        totalseg = totalseg+seg;
        totalhitSEG =totalhitSEG+1;
    end
    
    org = sum(sum(L==i));
    totalorg = totalorg+org;
end

precision = totalhit /(totalseg + (totalseg == 0));
recall = totalhit / (totalorg + (totalorg == 0));
penalty = totalhitSEG/totalhitGT;
