% merge all the small regions
C = bwconncomp(ucm < k);
labels = labelmatrix(C);
numRegions = length(unique(labels)) - 1;
S = regionprops(C); % region property in C
se1 = strel('square',3);
for i = 1: numRegions
    labels(imdilate(labels == i,se1)) = i;
end


ii = unique(labels);
%ii = ii(2:end); % remove zero label

colorsegs = ones([size(labels),size(im,3)]);
for i=1:length(ii)
    m = labels==ii(i);
    for c=1:size(im,3)
        tmp = imresize(im(:,:,c),1);
        colorsegs_tmp = colorsegs(:,:,c);
        colorsegs_tmp(m) = mean(tmp(m));
        colorsegs(:,:,c) = colorsegs_tmp;
        S(i).color(c) = mean(tmp(m));
    end
    hue = getFeatures(reshape(S(i).color,[1 1 3]),1,{'hue opp'},[]);
    S(i).hue = hue{1};
end

% find neighboring pixels
%d = 1:6; % offset
%glcm = zeros(numRegions+1,numRegions+1,8);
area_th = numel(labels)/20; % 1/20 of the total area
merge_flag = 1; % flage for keep merging
numMerges = 0;
while merge_flag
    merge_flag = 0;
    numMerges = numMerges + 1;
    fprintf('go through the loop %d times\n',numMerges);
    offset = [0 1; 0 -1; 1 0; -1 0];
    glcm = graycomatrix(labels,'Offset',offset, 'symmetric',true,'GrayLimits',[0 numRegions],...
        'NumLevels',numRegions+1);
    glcm = sum(glcm,3); % sum number of co-occurence
    neighbors = glcm > 0; % A(i,j) = 1 if i and j are neighbors
    % note that the first row and column are for 0 - the boundaries
    % so first row and first col are all 1, just take them out
    neighbors = neighbors(2:end,2:end);
    
    for i = 1: numRegions
        %include big ones to merging procedure, but omit really big guys
        %like >20000 , including them also works for some images
        if S(i).Area > area_th || S(i).Area == 0
            %         if S(i).Area == 0
            
            continue;
        end
        neighbor_i = find(neighbors(i,:) == 1);
        neighbor_i = neighbor_i(neighbor_i ~= i);
        if isempty(neighbor_i);
            fprintf('region %d does not have any neighbor?\n',i);
            continue;
        end;
        %     for j = 1:length(neighbor_i)
        %         %if neighbor_i(j) < i; continue; end
        %         if S(neighbor_i(j)).Area > area_th
        %            labels(labels == i) = neighbor_i(j);
        %            fprintf('Merge region %d into region %d\n',i,neighbor_i(j));
        %            break;
        %         end
        %     end
        [area,ind] = max(cat(1,S(neighbor_i).Area));
        %                 if (neighbor_i(ind) > i || area > area_th) && abs(circ_dist(S(i).hue,S(neighbor_i(ind)).hue)) < 1.25
        %look only for color similarity, do not look for other properties
        %threshold other than 1.25 might work as well but for that we need
        %to do parameter analysis, this 1.25 should be enough for now.
        if abs(circ_dist(S(i).hue,S(neighbor_i(ind)).hue)) < 1.25
            
            labels(labels == i) = neighbor_i(ind);
            %             fprintf('Merge region %d into region %d\n',i,neighbor_i(ind));
            S(neighbor_i(ind)).Area = S(neighbor_i(ind)).Area + S(i).Area;
            S(i).Area = 0;
            merge_flag = 1;
        end
    end
    
    
end

%merge really small ones (<20) to it's largest neighbor
for i = 1: numRegions
    %         if S(i).Area > area_th || S(i).Area == 0
    if S(i).Area == 0
        continue;
    end
    neighbor_i = find(neighbors(i,:) == 1);
    neighbor_i = neighbor_i(neighbor_i ~= i);
    if isempty(neighbor_i);
        continue;
    end;
    %     for j = 1:length(neighbor_i)
    %         %if neighbor_i(j) < i; continue; end
    %         if S(neighbor_i(j)).Area > area_th
    %            labels(labels == i) = neighbor_i(j);
    %            fprintf('Merge region %d into region %d\n',i,neighbor_i(j));
    %            break;
    %         end
    %     end
    [area,ind] = max(cat(1,S(neighbor_i).Area));
    %         if (neighbor_i(ind) > i || area > area_th) && abs(circ_dist(S(i).hue,S(neighbor_i(ind)).hue)) < 1.25
    if S(i).Area < 20
        
        labels(labels == i) = neighbor_i(ind);
        %             fprintf('Merge region %d into region %d\n',i,neighbor_i(ind));
        S(neighbor_i(ind)).Area = S(neighbor_i(ind)).Area + S(i).Area;
        S(i).Area = 0;
        merge_flag = 1;
    end
end