%

function [A,mdist]= fast_buildW(im,d_max, p, rf, mixture_params, opts)

which_feature = opts.features.which_features;
which_affinity = opts.affinityFunction;
%minAffty=.01;	% Impose a lower bound on affinity! helps
		% remove single pixel regions.

sizeIm=size(im);
A=speye(sizeIm(1)*sizeIm(2));

im=reshape(im,sizeIm(1)*sizeIm(2),1);

for i=0:d_max
 for j=-d_max:d_max

  if (j>0 || i>0)
   % Figure out which diagonal we are working on
   dn=(sizeIm(1)*i)+j;
   fprintf(2,'j=%d, i=%d, dn=%d\n',j,i,dn);
 
   % left side of the main diagonal
   % Vector for pixels along the main diagonal
   indx1 = 1:length(im)-dn; % left side
   indx2 = 1+dn:length(im); % right side
   Mask = ones(length(indx1),1);
   % build the mask where the values are filled in
   if (j>0)
    for k=1:length(Mask)
     if (mod(k,sizeIm(1))==0||mod(k,sizeIm(1))>sizeIm(1)-j); Mask(k)=0; end;
    end;    
   else if (j<0)
     for k=1:length(Mask)
      if (mod(k,sizeIm(1))>0 && mod(k,sizeIm(1))<=abs(j)); Mask(k)=0; end;
     end;
    end;
   end;   
   
   if strcmp(which_affinity,'difference')
       if strcmp(which_feature,'luminance')
           F1 = im(indx1,:); 
           F2 = im(indx2,:);
           Fdist = ((F1-F2) .^2).*Mask;    
       end
   end
   %%
   diag1=[zeros((sizeIm(1)*sizeIm(2))-length(Fdist),1)' Fdist']';
   diag2=[Fdist' zeros((sizeIm(1)*sizeIm(2))-length(Fdist),1)']';

   % Stuff the appropriate diagonals in A
   A=spdiags(diag1,dn,A);
   A=spdiags(diag2,-dn,A);
  end;

 end;
end;

rho     = 1.5;
  
mdist = median(sqrt(nonzeros(A)));
sigma = rho*mdist; 
A = exp(-A./(2 * sigma * sigma));



%im=reshape(im,sizeIm(1),sizeIm(2));

