
function [precision,recall] = evalPrecisionRecall(groundTruth,result)

gto = groundTruth{1,1};

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
if ~isempty(whiteLabel)
    noback = gt~=stromaLabel &  gt~=whiteLabel;
elseif ~isempty(stromaLabel)
    noback =  gt~=stromaLabel;
else
    noback = ones(size(gt));
end
L= bwlabeln(noback);

totalRegion=zeros(max(max(result)),max(max(L)));
for i = 1:size(result,1)
    for j = 1:size(result,2)
        rno= result(i,j);
        gs = L(i,j);
        if gs~=0 && rno~=1
            totalRegion(rno,gs)=totalRegion(rno,gs)+1;
        end
    end
end

totalhit =0;
totalseg=0;
totalorg = 0;
for i=1:max(max(L))
    [hit,in] = max(totalRegion(:,i));
    totalhit = totalhit +hit;
    seg = sum(sum(result==in));
    totalseg = totalseg+seg;
    org = sum(sum(L==i));
    totalorg = totalorg+org;
end

precision = totalhit / totalseg;
recall = totalhit / totalorg;