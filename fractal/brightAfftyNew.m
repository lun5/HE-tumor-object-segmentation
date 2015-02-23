% What happens if we get affinities in a larger neighborhood?

function [A]=brightAfftyNew(im,d_max,mDist)

minAffty=.01;	% Impose a lower bound on affinity! helps
		% remove single pixel regions.

sizeIm=size(im);
maxNZ=sizeIm(1)*sizeIm(2)*(((d_max*2)+1)^2);
A=speye(sizeIm(1)*sizeIm(2));

im=reshape(im,sizeIm(1)*sizeIm(2),1);

for i=0:d_max
 for j=-d_max:d_max

  if (j>0 || i>0)
   % Figure out which diagonal we are working on
   dn=(sizeIm(1)*i)+j;
   %fprintf(2,'j=%d, i=%d, dn=%d\n',j,i,dn);
 
   % left side of the main diagonal
   % Vector for pixels along the main diagonal
%   F1=[im(1:length(im)-dn) j*ones(length([1+dn:length(im)]),1) i*ones(length([1+dn:length(im)]),1)];
   F1=[im(1:length(im)-dn) zeros(length([1+dn:length(im)]),1) zeros(length([1+dn:length(im)]),1)];
   F1 = double(F1); % Luong added 
   % Corresponding pixels vector
   F2=[im(1+dn:length(im)) zeros(length([1+dn:length(im)]),1) zeros(length([1+dn:length(im)]),1)];
   F2 = double(F2); % Luong added
   Mask=ones(length(F2),1);

   if (j>0)
    for k=1:length(Mask)
     if (mod(k,sizeIm(1))==0||mod(k,sizeIm(1))>sizeIm(1)-j) Mask(k)=0; end;
    end;    
   else if (j<0)
     for k=1:length(Mask)
      if (mod(k,sizeIm(1))>0 && mod(k,sizeIm(1))<=abs(j)) Mask(k)=0; end;
     end;
    end;
   end;   
   
   %Fdist=sqrt(sum((mod(F1-F2,255)).^2,2)); 
   Fdist=sqrt(sum((F1-F2).^2,2)); 
   %Fdist = abs(sum(F1-F2,2));
   Fdist=(exp(-(Fdist.^2)/(2*mDist^2))+minAffty).*Mask;
%   Fdist=(exp(-(Fdist.^2)/(2*mDist^2))+eps).*Mask;
%   Fdist=exp(-(Fdist)/(2*mDist^2))+eps.*Mask;
   diag1=[zeros((sizeIm(1)*sizeIm(2))-length(Fdist),1)' Fdist']';
   diag2=[Fdist' zeros((sizeIm(1)*sizeIm(2))-length(Fdist),1)']';

   % Stuff the appropriate diagonals in A
   A=spdiags(diag1,dn,A);
   A=spdiags(diag2,-dn,A);
  end;

 end;
end;


%im=reshape(im,sizeIm(1),sizeIm(2));

