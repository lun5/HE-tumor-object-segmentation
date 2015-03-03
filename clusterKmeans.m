function [discComp,U,S] = clusterKmeans(A,dm,U,S)
%function [discComp,U,S] = clusterKmeans(Pts,A,dm,U,S)

  FALSE = (0 == 1);
  TRUE = ~FALSE;

  Pts = A;
  
  % variable dm gets overwritten at the end
  dmo = dm;
  
  if (nargin == 2) % or 3
    
    %A = A/median((sum(A,1)));
    
    D = sum(A, 1)';              % Normalize column sum to one. 
    sqrtDinv = spdiags(D.^(-0.5),0,length(D),length(D));%(sqrtD .^ -1) * ones(1, length(D)); 
    Mcut = sqrtDinv * A * sqrtDinv;         % M = D^-0.5 Markov D^0.5 

    [U,S,V] = svds(Mcut,dm);
    %[U S V] = svd(Mcut); 
    S = diag(S);

    %[U,S] = eigs(Mcut,20,'LM');    
    
    clear Q
    clear V
  end
  %keyboard;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %%%%  Projection Clustering Using K-means %%%%%%%%%%%%%%%%%%%%% 
  %%%%  See Ng, Jordan, and Wiess           %%%%%%%%%%%%%%%%%%%%% 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  %% Params for selecting initial K-means centers 
  minRes2 = 0.25;   
  MARKOV = FALSE; % MARKOV = TRUE; minRes2 = 0.5; nIts = 64; 
  
  
  %%Rerun: rand('state', saveRandState); 
  
  n = size(Pts); 
  nIts  = 2;
  %% Set up data for center selection: 
  if MARKOV 
    c = U .* (ones(n(1), 1) * (S .^ nIts)'); 
    E = c; 
    E = ((sum(E.*E, 2).^-0.5) * ones(1, size(E,2))) .* E; 
  else 
    E = U(:, 1:dm); 
    E = ((sum(E.*E, 2).^-0.5) * ones(1, dm)) .* E; 
  end 
  
  %% D0 center selection: 
  %% Rerun: rand('state', saveRandState); 
  saveRandState = rand('state'); 
  
  %% Select initial centers to be nearly orthogonal. 
  j = round(0.5 + n(1) * rand(1)); 
  j = max(j, 1); j = min(j,n(1)); 
  c = E(j,:); 
  res = E; 
  k = 2; 
  while MARKOV | k <= dm 
    res = res - (res * (c(k-1,:)')) * c(k-1,:); 
    nrmRes2 = sum(res .* res, 2); 
    samp = cumsum(nrmRes2 .* (nrmRes2>minRes2)); 
    if samp(n(1)) == 0 | isnan(samp(n(1))) | isinf(samp(n(1))) 
      break; 
    end  
    samp = samp/samp(n(1)); 
    r = rand(1); 
    idx = find(samp>=r); 
    if any(idx) 
      c = [c ; E(idx(1), :)]; 
      k = k+1; 
    else 
      error('Random draw fanned!??!'); 
    end  
  end 
  k = k-1; 
  if k < dm & ~MARKOV 
    fprintf(2,' Got only %d basis elements\n', k); 
    dm = k; 
  else 
    if (k ~= dmo)
      fprintf(2,' Got %d basis elements\n', k); 
    end
    dm = k; 
  end 
  
  %% Call kmeans 
  options = foptions; 
  [centers options post errlog] = kmeanslocal(c, E, options); 
  discComp = post>0; 
