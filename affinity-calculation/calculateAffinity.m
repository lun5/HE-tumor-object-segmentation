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
function [Pts,A,mdist] = calculateAffinity(I,opts)
        
    %% calculate features
    which_features = opts.features.which_features;
    f_maps = [];
    for feature_iter = 1: length(which_features)
        if (opts.display_progress), fprintf('\nProcessing feature type ''%s'':\n',which_features{feature_iter}); end
        f_maps = cat(3,f_maps,getFeatures(double(I),1,which_features{feature_iter},opts));
    end
    
    %% calculate affinity
    rho     = 1.5;%afftyPar.rho;
    sizeIm  = size(f_maps);%afftyPar.sizeIm;
    if (opts.display_progress), fprintf('\nProcessing affinity function ''%s'':\n',opts.affinityFunction); end
    
    d_max = opts.localPairs.rad;
    A=speye(sizeIm(1)*sizeIm(2));
    Mask_mt = speye(sizeIm(1)*sizeIm(2));
    
    %% Type of features and affinity functions
    which_affinity = opts.affinityFunction;
    
    %% learn joint density in case of PMI
    if strcmp(which_affinity,'PMI')
     for feature_iter = 1: length(which_features)
      if strcmp(which_features(feature_iter),'luminance')        
        [~, index] = ismember('luminance', which_features);
        p_luminance = learnP_A_B(f_maps(:,:,index),opts);
        %% learn w predictor
        if (opts.approximate_PMI)
            rf_luminance = learnPMIPredictor(f_maps(:,:,index),p_luminance,opts);
        end
      elseif strcmp(which_features(feature_iter),'brightness opp')
        [~, index] = ismember('brightness opp', which_features);
        p_bright = learnP_A_B(f_maps(:,:,index),opts);
        %% learn w predictor
        if (opts.approximate_PMI)
            rf_bright = learnPMIPredictor(f_maps(:,:,index),p_bright,opts);
        end
      elseif strcmp(which_features(feature_iter),'saturation opp')
        [~, index] = ismember('saturation opp', which_features);
        p_sat = learnP_A_B(f_maps(:,:,index),opts);
        %% learn w predictor
        if (opts.approximate_PMI)
            rf_sat = learnPMIPredictor(f_maps(:,:,index),p_sat,opts);
        end
      elseif strcmp(which_features(feature_iter),'hue opp');
        [~, index] = ismember('hue opp', which_features);
        Nsamples = 10000;
        F = sampleF(f_maps(:,:,index),Nsamples,opts);  
        [ params,~, prior_probs] = mixture_of_bivariate_VM(F, 6);
        mixture_params.params = params;
        mixture_params.prior_probs = prior_probs;
      end
     end
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
   indx1 = 1:sizeIm(1)*sizeIm(2) - dn;%1:length(im)-dn;
   indx2 = 1 + dn:sizeIm(1)*sizeIm(2);%1+dn:length(im);
   Fdist = ones(length(indx1),1); 
   for feature_iter = 1:length(which_features)
     im=reshape(f_maps(:,:,feature_iter),sizeIm(1)*sizeIm(2),1);
     F1=[im(indx1)];% zeros(length([1+dn:length(im)]),1) zeros(length([1+dn:length(im)]),1)];
     F1 = double(F1); % Luong added 
     % Corresponding pixels vector
     F2=[im(indx2)];% zeros(length([1+dn:length(im)]),1) zeros(length([1+dn:length(im)]),1)];
     F2 = double(F2); % Luong added
     % where it went wrong (16,1), (32,17)

     %Fdist=sqrt(sum((mod(F1-F2,255)).^2,2)); 
     %% different type of features and affinity
     if strcmp(which_affinity,'difference')
      if ismember(which_features,'luminance')
        Fdist = sum((F1-F2).^2,2);
      elseif ismember(which_features,'hue opp')
        Fdist = circ_dist(F1(:,1),F2(:,1));
      end
   % calculate PMI in case of PMI. What do I do with difference?
     elseif strcmp(which_affinity,'PMI')
       if strcmp(which_features{feature_iter},'luminance')
        if (opts.approximate_PMI)
            pmi = fastRFreg_predict([F1 F2],rf_luminance);
        else
            [pmi,~,~] = evalPMI(p_luminance,[F1 F2],[],[],[],opts);            
        end
       elseif strcmp(which_features{feature_iter},'brightness opp')
        if (opts.approximate_PMI)
            pmi = fastRFreg_predict([F1 F2],rf_bright);
        else
            [pmi,~,~] = evalPMI(p_bright,[F1 F2],[],[],[],opts);
        end
       elseif strcmp(which_features{feature_iter},'saturation opp')
        if (opts.approximate_PMI)
            pmi = fastRFreg_predict([F1 F2],rf_sat);
        else
            [pmi,~,~] = evalPMI(p_sat,[F1 F2],[],[],[],opts);
        end
       elseif strcmp(which_features{feature_iter},'hue opp')
        [pmi,~,~] = evalPMI_theta([F1 F2], mixture_params, opts); 
        %pmi = log(pmi);
       end
       Fdist = Fdist.*pmi;
    end
   end
   Mask=ones(length(F2),1);
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

   Fdist=Fdist.*Mask; 
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
Mask_mt = Mask_mt - speye(sizeIm(1)*sizeIm(2));
[row,col] = find(Mask_mt);
lin_indx = sub2ind(size(Mask_mt),row,col);
data = full(A(lin_indx));
clear A Mask_mt;

if strcmp(which_affinity,'difference')
  if strcmp(which_features,'luminance')
    mdist = median(sqrt(nonzeros(data)));
    %Fdist = abs(sum(F1-F2,2));
    sigma = rho*mdist;
    aff = exp(-data/(2*sigma^2));
    A = sparse(row,col,aff,prod(sizeIm),prod(sizeIm))+ speye(sizeIm(1)*sizeIm(2));
    Pts(:,1) = f_maps(:).*255;
    clear data aff;
  elseif strcmp(which_features,'hue opp')
    thetahat = circ_mean(nonzeros(data)); mdist = circ_kappa(nonzeros(data));
    Pts(:,1) = (f_maps(:) + pi)/(2*pi)*255;  
    aff = exp(mdist*cos(data-thetahat))./exp(mdist);
    A = sparse(row,col,aff,prod(sizeIm),prod(sizeIm))+ speye(sizeIm(1)*sizeIm(2));
    clear data aff;
  end
elseif strcmp(which_affinity,'PMI')
  % normalize to between 0 and 1 - DO WE NEED THIS?
  prc = 2; lb = prctile(data,prc); ub = prctile(data,100- prc);
  aff = data; aff(aff<lb) = lb;
  aff = (aff -lb)./(ub - lb);
  %aff = (data - min(data))./(max(data) - min(data));
  A = sparse(row,col,aff,sizeIm(1)*sizeIm(2),sizeIm(1)*sizeIm(2))+ max(aff)*speye(sizeIm(1)*sizeIm(2));
  f_maps_curr = f_maps(:,:,1);
  Pts(:,1)= (f_maps_curr(:) - min(f_maps_curr(:)))./(max(f_maps_curr(:))-min(f_maps_curr(:)))*255;
  mdist = [];
end
      
end