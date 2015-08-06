%% map label to output color to use in TextonBoost program
% Luong Nguyen 8/5/2015
% Input: label (id from 1-> whatever) where each id is associated with one
% label: e.g. "carcinoma", etc
% Output: col = [rgb] value using similar code with the C# ==> such bs
function col = mapLabel2Color(label)
    r = 0; g = 0; b = 0; label = uint8(label);
    for j = 0:7
       r = bitor(r,bitshift(bitand(bitsra(label,0),1),7-j));
       g = bitor(g,bitshift(bitand(bitsra(label,1),1),7-j));
       b = bitor(b,bitshift(bitand(bitsra(label,2),1),7-j));
       label = bitsra(label,3);
       if label <= 0
           break;
       end
    end
    col = [r g b];
end