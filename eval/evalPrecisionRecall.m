
function [precision,recall, penalty, ri, gce, vi] = evalPrecisionRecall(groundTruth,result)

if max(result(:)) <= 1
    precision = 0;
    recall = 0;
    penalty = 0;
    ri=0;
    vi=10;
    gce=1;
    return;
end

gto = groundTruth{1,1};

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
logical_cells = cellfun(cellfind('stroma'),gto.names);
len= length(logical_cells)+1;
stromaLabel=len;
if sum(logical_cells) > 0
    stromaLabel = find(logical_cells==1);
end
whiteLabel = len;
logical_cells = cellfun(cellfind('white'),gto.names);
if sum(logical_cells) > 0
    whiteLabel = find(logical_cells==1);
end
gt = gto.Segmentation;
gtorig = gt;
gt(gt==0)=stromaLabel;
noback = gt~=stromaLabel &  gt~=whiteLabel;
noback(gto.Boundaries)=0;
L= bwlabeln(noback);

totalRegion=zeros(max(max(result)),max(max(L))+1);
for i = 1:size(result,1)
    for j = 1:size(result,2)
        rno= result(i,j);
        gs = L(i,j);
        % if gs~=0
        totalRegion(rno,gs+1)=totalRegion(rno,gs+1)+1;
        %end
    end
end

totalhit =0;
totalseg=0;
totalorg = sum(sum(L>0));
totalhitGT = 0;
totalhitSEG=0;
totalhitSEGarea=0;
processedGT =zeros(1,max(max(L))+1);
processedGTArea = zeros(1,max(max(L))+1);
countSEGperGT= zeros(1,max(max(L))+1);
countHitSeg=0;
hitGT =zeros(1,max(max(L))+1);
for i=1:size(totalRegion,1)
    flag=true;
    segFlag=false;
    prevMax=0;
    while flag
        indMax = find(totalRegion(i,:)== max(totalRegion(i,:)));
        
        if indMax==1
            TT=totalRegion(i,indMax)./sum(totalRegion(i,:));
            if TT>.88
                segFlag=true;
                break;
            else
                totalRegion(i,indMax)=-1;
                continue;
            end
        end
        
        if sum(totalRegion(i,indMax))<=0
            if prevMax >10
                totalhitSEG =totalhitSEG+1;
                totalhitSEGarea = totalhitSEGarea + prevMax;
                countHitSeg=countHitSeg+1;
                segFlag=true;
                break;
            else
                segFlag=true;
                break;
            end
        else
            indMax = indMax(1);
        end
        
        
        prevMax = max(totalRegion(i,:));
        if processedGT(1,indMax)==1
            totalRegion(i,indMax)=-1;
        else
            flag=false;
        end
        
    end
    if segFlag
        continue;
    end
    if totalRegion(i,indMax) < 10
        continue;
    end
    
    %     GToriglabel = round(mean(mean(gtorig (result==i))));
    %     if GToriglabel==stromaLabel
    %         continue;
    %     elseif GToriglabel==whiteLabel
    %         continue;
    %     end
    
    
    
    countHitSeg=countHitSeg+1;
    hit = totalRegion(i,indMax);
    totalhit = totalhit +hit;
    seg = sum(sum(result==i));
    totalseg = totalseg+seg;
    totalhitSEGarea = totalhitSEGarea + seg;
    
    if hitGT(1,indMax)==0
        hitGT(1,indMax)=1;
        totalhitGT = totalhitGT+1;
        
    end
    
    thr = sum(sum(L==indMax-1))/10;%20; %CHANGE TO 10 to see if can save penalty
    if hit > sum(sum(L==indMax-1)) - thr
        processedGT(1,indMax)=1;
        countSEGperGT(indMax)=countSEGperGT(indMax)+1;
    else
        processedGTArea(indMax)= processedGTArea(indMax) +hit;
        countSEGperGT(indMax)=countSEGperGT(indMax)+1;
        if processedGTArea(indMax) > sum(sum(L==indMax-1)) - thr
            processedGT(1,indMax)=1;
        end
    end
end

if totalseg==0
    precision=0;
    penalty=1;
else
    precision = totalhit / totalseg;
    penalty = sqrt(sum(countSEGperGT))/sqrt(totalhitGT); %
end

if totalorg ==0
    recall=0;
else
    recall = totalhit / totalorg;
end

[ri,gce,vi]=compare_segmentations(result,L);


