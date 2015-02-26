%

function [A]= fast_buildW(im,d_max, p, rf, mixture_params, opts)

which_feature = opts.features.which_features;
which_affinity = opts.affinityFunction;
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
   fprintf(2,'j=%d, i=%d, dn=%d\n',j,i,dn);
 
   % left side of the main diagonal
   % Vector for pixels along the main diagonal
%   F1=[im(1:length(im)-dn) j*ones(length([1+dn:length(im)]),1) i*ones(length([1+dn:length(im)]),1)];
   indx1 = 1:length(im)-dn;
   indx2 = 1+dn:length(im);
   %F1= [im(1:length(im)-dn)]; %zeros(length([1+dn:length(im)]),1) zeros(length([1+dn:length(im)]),1)];
   %F1 = double(F1); % Luong added 
   % Corresponding pixels vector
   %F2= [im(1+dn:length(im))]; % zeros(length([1+dn:length(im)]),1) zeros(length([1+dn:length(im)]),1)];
   %F2 = double(F2); % Luong added
   %Mask=ones(length(F2),1);
   Mask = ones(length(indx1),1);
   
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
       % make ht number of layers in pyramid
       ht = opts.pyramid_ht;
       [pyr, pind] = buildLpyr(im,ht);
       im_pyr = zeros(sizeIm(1)*sizeIm(2),ht);
       for pyr_iter = 1:ht
           res = reconLpyr(pyr, pind, pyr_iter);
           im_pyr(:,pyr_iter) = res;%reshape(res,sizeIm(1)*sizeIm(2),1);
       end
       F1 = im_pyr(indx1,:); 
       F2 = im_pyr(indx2,:);
       
        %Fdist=sqrt(sum((mod(F1-F2,255)).^2,2)); 
        %Fdist=sqrt(sum(mod(F1-F2,pi).^2,2)); 
        Fdist = mod(abs(F1-F2),pi);
        %Fdist = abs(sum(F1-F2,2));
        mDist = median(Fdist);
        Fdist=(exp(-(Fdist.^2)/(1.5^2*mDist^2))+minAffty).*Mask;
        %   Fdist=(exp(-(Fdist.^2)/(2*mDist^2))+eps).*Mask;
        %   Fdist=exp(-(Fdist)/(2*mDist^2))+eps.*Mask;
   elseif strcmp(which_affinity,'PMI')
       F1 = im(indx1); F2 = im(indx2);
       if strcmp(which_feature,'luminance')
           if isempty(p) && isempty(rf)
               error('prog:input',...
                   'missing input values for affinity calculation with feature %s and affinity function %s',...
                   which_feature, which_affinity);
           end
           if (opts.approximate_PMI)
               Fdist = fastRFreg_predict([F1 F2],rf);
           else
               [Fdist,~,~] = evalPMI(p,[F1 F2],[],[],[],opts);
           end
           %reg = prctile(nonzeros(Fdist),5);
           %Fdist = exp(Fdist);%+reg);
           Fdist = Fdist - min(Fdist) + 0.01;
           Fdist = Fdist.*Mask;
       elseif strcmp(which_feature,'hue opp')
           [Fdist,~,~] = evalPMI_theta([F1 F2], mixture_params, opts);
           Fdist = Fdist.*Mask;
           %reg = prctile(nonzeros(Fdist),5);
           %Fdist = log(evalPMI_theta+reg);
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


%im=reshape(im,sizeIm(1),sizeIm(2));

