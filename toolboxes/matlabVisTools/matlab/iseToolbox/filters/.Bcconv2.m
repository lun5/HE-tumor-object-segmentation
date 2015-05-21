function result = cconv2(im,filter)
% CCONV2: Circular convolution (calls upConv).
%
% res = cconv2(image,filter)
%
% DJH, 8/96
% update 12/97

result = upConv(im,filter,'circular',[1,1]);
return;

%%% Debug

%create a random noise image
im=round(rand(100,100));

%create a filter.
filter = [1,2,4,2,1]'*[1,2,4,2,1];
filter=filter/sum(sum(filter));

res=cconv2(im,filter);
