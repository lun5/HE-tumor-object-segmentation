%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%  Projection Clustering Using K-means %%%%%%%%%%%%%%%%%%%%% 
%%%%  See Ng, Jordan, and Wiess           %%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% U has eigenvectors
%% S has eigenvalues
 
%% Params for selecting initial K-means centers 
minRes2 = 0.25;   
MARKOV = false; % MARKOV = TRUE; minRes2 = 0.5; nIts = 64; 
 
 
%%Rerun: rand('state', saveRandState); 
 
%n = size(Pts); 
n = size(U,1);
%% Set up data for center selection: 

% set up the k-dim unit sphere
E = U(:, 1:dm); 
E = ((sum(E.*E, 2).^-0.5) * ones(1, dm)) .* E; 
 
 
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
  fprintf(2,' Got %d basis elements\n', k); 
  dm = k; 
end 
 
%% Call kmeans
options = foptions; 
[centers options post errlog] = kmeans(c, E, options); 
discComp = post>0; 
