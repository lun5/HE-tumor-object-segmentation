function m=initializeClusterVectorsUsingPrincipalComponent(Data , k )
    CR = corrcoef(Data);
    [V,D] = eig(CR);
    newdata = Data * V;
    maxentry = max(newdata(:,1));
    minentry = min(newdata(:,1));
    d = (maxentry - minentry)/k;
    c = zeros(k+1);
    c(1) = minentry;
    for i = 2 : k+1
        c(i) = c(i - 1) + d;
    end
    No = zeros(1,k);
    newM = zeros(k, size(Data,2));
    [rows, cols] = size(newdata);
    for i = 1 : rows
        for j = 1 : k 
            if ((newdata(i, 1) >= c(j)) && (newdata(i, 1) < c(j + 1))) 
                for l = 1 :cols
                    newM(j, l)= newM(j, l) + newdata(i, l);
                end
                No(j) = No(j) + 1;
                break;
            end
        end
    end
    for i = 1 : k
        for j = 1 : cols
            newM(i, j)= newM(i, j) / No(i);
        end
    end
%     temp = inv(V);
    m = newM/V;
end