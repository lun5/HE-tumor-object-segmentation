%% map output color to label to evaluate TextonBoost performance
% Luong Nguyen 8/6/2015
% Input: col = [rgb] value using similar code with the C# ==> such bs
% Output: label (id from 1-> whatever) where each id is associated with one
% label: e.g. "carcinoma", etc

function label = mapColor2Label(col)
    r = col(:,1); g = col(:,2); b = col(:,3); 
    label = uint8(0);
    for j = 0:7
        label = bitor(bitor(bitshift(label,3),bitshift(bitand(bitsra(r,j),1),0)),...
            bitor(bitshift(bitand(bitsra(g,j),1),1),bitshift(bitand(bitsra(b,j),1),2)));
    end
end