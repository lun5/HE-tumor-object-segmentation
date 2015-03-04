%% function [Ws,im_sizes] = getW(I,opts_affinity)
% builds affinity matrices Ws for image I
% 
% INPUTS
%  I                - NxMxC query image
%  opts_affinity    - parameter settings (see setEnvironment_affinity)
%
% OUTPUTS
%  Ws         - affinity matrices; Ws{i} is the affinity matrix for the image at scale i
%  im_sizes   - im_sizes{i} gives the dimensions of the image at scale i
%               (note: dimensions are num cols x num rows; this is the
%                opposite of matlab's default!)
%
% -------------------------------------------------------------------------
% HE segmentation Toolbox
% Luong Nguyen, 2014 [lun5@pitt.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------
% taken from 
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2014 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------
function [Pts,A,mdist] = calculateAffinity(I,opts)
        
    %% calculate features
    if (opts.display_progress), fprintf('\nProcessing feature type ''%s'':\n',opts.features.which_features{1}); end
    f_maps = getFeatures(double(I),1,opts.features.which_features,opts);
    
    %% calculate affinity
    rho     = 1.5;%afftyPar.rho;
    sizeIm  = size(f_maps);%afftyPar.sizeIm;
    if (opts.display_progress), fprintf('\nProcessing affinity function ''%s'':\n',opts.affinityFunction); end
    im=reshape(f_maps,sizeIm(1)*sizeIm(2),1);
    d_max = opts.localPairs.rad;
    A=speye(sizeIm(1)*sizeIm(2));
    Mask_mt = speye(sizeIm(1)*sizeIm(2));
    
    %% Type of features and affinity functions
    which_features = opts.features.which_features{1};
    which_affinity = opts.affinityFunction;
    
    if strcmp(which_features,'luminance') && strcmp(which_affinity,'PMI')
        p = learnP_A_B(f_maps,opts);
        %% learn w predictor
        if (opts.approximate_PMI)
            rf = learnPMIPredictor(f_maps_curr,p,opts);
        end
    elseif strcmp(which_features,'hue opp') && strcmp(which_affinity,'PMI');
        Nsamples = 10000;
        F = sampleF(f_maps_curr,Nsamples,opts);  
        [ params,~, prior_probs] = mixture_of_bivariate_VM(F, 6);
        mixture_params.params = params;
        mixture_params.prior_probs = prior_probs;
    end
 
for i=0:d_max
 for j=-d_max:d_max

   if (j>0 || i>0)
   % Figure out which diagonal we are working on
   dn=(sizeIm(1)*i)+j;
   %fprintf(2,'j=%d, i=%d, dn=%d\n',j,i,dn);
 
   % left side of the main diagonal
   % Vector for pixels along the main diagonal
   %F1=[im(1:length(im)-dn) j*ones(length([1+dn:length(im)]),1) i*ones(length([1+dn:length(im)]),1)];
   F1=[im(1:length(im)-dn)];% zeros(length([1+dn:length(im)]),1) zeros(length([1+dn:length(im)]),1)];
   F1 = double(F1); % Luong added 
   % Corresponding pixels vector
   F2=[im(1+dn:length(im))];% zeros(length([1+dn:length(im)]),1) zeros(length([1+dn:length(im)]),1)];
   F2 = double(F2); % Luong added
   Mask=ones(length(F2),1);
   % where it went wrong (16,1), (32,17)
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
   %% different type of features and affinity
   if strcmp(which_affinity,'difference')
      if strcmp(opts.features.which_features{1},'luminance')
        Fdist = sum((F1-F2).^2,2);
      elseif strcmp(opts.features.which_features{1},'hue opp')
        Fdist = circ_dist(F1(:,1),F2(:,1));
      end
   elseif strcmp(which_affinity,'PMI')
      if strcmp(opts.features.which_features{1},'luminance')
        if (opts.approximate_PMI)
            Fdist = fastRFreg_predict([F1 F2],rf);
        else
            [Fdist,~,~] = evalPMI(p,[F1 F2],[],[],[],opts);
        end 
      elseif strcmp(opts.features.which_features{1},'hue opp')
        [Fdist,~,~] = evalPMI_theta([F1 F2], mixture_params, opts);
      end
   end
   
   %Fdist(Fdist == 0) = minAffty;
   Fdist=Fdist.*Mask; 
   %Fdist=(exp(-(Fdist.^2)/(2*mDist^2+minAffty))).*Mask;
%   Fdist=(exp(-(Fdist.^2)/(2*mDist^2))+eps).*Mask;
%   Fdist=exp(-(Fdist)/(2*mDist^2))+eps.*Mask;
   diag1=[zeros((sizeIm(1)*sizeIm(2))-length(Fdist),1)' Fdist']';
   diag2=[Fdist' zeros((sizeIm(1)*sizeIm(2))-length(Fdist),1)']';

   % Stuff the appropriate diagonals in A
   A=spdiags(diag1,dn,A);
   A=spdiags(diag2,-dn,A);
   
   diag1_m=[zeros((sizeIm(1)*sizeIm(2))-length(Mask),1)' Mask']';
   diag2_m=[Mask' zeros((sizeIm(1)*sizeIm(2))-length(Mask),1)']';
   Mask_mt = spdiags(diag1_m,dn,Mask_mt);
   Mask_mt=spdiags(diag2_m,-dn,Mask_mt);
  end;

 end;
end;
A = A - speye(sizeIm(1)*sizeIm(2));
[row,col] = find(Mask_mt);
lin_indx = sub2ind(size(Mask_mt),row,col);
data = full(A(lin_indx));clear A Mask_mt;

if strcmp(which_affinity,'difference')
  if strcmp(which_features,'luminance')
    mdist = median(sqrt(nonzeros(data)));
    %Fdist = abs(sum(F1-F2,2));
    sigma = rho*mdist;
    aff = exp(-data/(2*sigma^2));
    A = sparse(row,col,aff,prod(sizeIm),prod(sizeIm));
    Pts(:,1) = f_maps(:).*255;
    clear data aff;
  elseif strcmp(which_features,'hue opp')
    thetahat = circ_mean(nonzeros(data)); mdist = circ_kappa(nonzeros(data));
    Pts(:,1) = (f_maps(:) + pi)/(2*pi)*255;  
    aff = exp(mdist*cos(data-thetahat))./exp(mdist);
    A = sparse(row,col,aff,prod(sizeIm),prod(sizeIm));
    clear data aff;
  end
elseif strcmp(which_affinity,'PMI')
  % normalize to between 0 and 1 - DO WE NEED THIS?
  aff = (data - min(data))./(max(data) - min(data));
  A = sparse(row,col,aff,prod(sizeIm),prod(sizeIm));
  Pts(:,1)= (f_maps(:) - min(f_maps(:)))./(max(f_maps(:))-min(f_maps(:)))*255;
  mdist = [];
end


%     dsThres =  opts.localPairs.rad + .1;%afftyPar.dsThres;
%     rho     = 1.5;%afftyPar.rho;
%     sizeIm  = size(f_maps);%afftyPar.sizeIm;
%     
%     Pts(:,1) = f_maps(:);
%     n = size(Pts);
%     [xi, yi] = meshgrid(1:sizeIm(2), 1:sizeIm(1));
%     Pts(:,2) = xi(:);
%     Pts(:,3) = yi(:);
%   
%     %% Distance between consecutive points, and median distance
%     X = Pts(:,1) * ones(1, n(1));
%     X = X - X';
%     Y = Pts(:,2) * ones(1, n(1));
%     Y = Y - Y';
%     Z = Pts(:,3) * ones(1, n(1));
%     Z = Z - Z';
%   
%     dS  = max(abs(Y),abs(Z));
%     d = (dS < dsThres).* (X .^2);
%     idx = find(d > 0);
%     mdist = median(sqrt(d(idx)));
%     sigma = rho*mdist; 
%     binNhbr = (dS < dsThres) & (dS ~= 0);
%     A = (dS < dsThres).*exp( -d/(2 * sigma * sigma));
      
end