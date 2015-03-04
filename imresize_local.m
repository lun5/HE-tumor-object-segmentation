function reconsIm = imresize_local(raw_image,lvl)
r = double(raw_image(:,:,1))/255; g = double(raw_image(:,:,2))/255; 
b = double(raw_image(:,:,3))/255; 

[pyrR,pindR] = buildGpyr(r, lvl);
[pyrG,pindG] = buildGpyr(g, lvl);
[pyrB,pindB] = buildGpyr(b, lvl);
reconsR = pyrR(end-prod(pindR(lvl,:))+1:end);
reconsR = reshape(reconsR,pindR(lvl,:))./2^(lvl-1);
reconsG = pyrG(end-prod(pindG(lvl,:))+1:end);
reconsG = reshape(reconsG,pindG(lvl,:))./2^(lvl-1);
reconsB = pyrB(end-prod(pindB(lvl,:))+1:end);
reconsB = reshape(reconsB,pindB(lvl,:))./2^(lvl-1);
reconsIm = cat(3,reconsR,reconsG,reconsB);
end