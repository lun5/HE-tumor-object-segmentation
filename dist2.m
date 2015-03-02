function d2 = dist2(data, centres)
numdata = size(data,1); 
numcenters = size(centres,1);
d2 = zeros(numdata, numcenters);
for i = 1:numdata
    for j = 1:numcenters
        d2(i,j) = sqrt(sum((data(i,:) - centres(j,:)).^2,2));
    end
end
end