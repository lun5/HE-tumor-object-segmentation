function [Pts,A,binNhbr,binSupp, mdist] = mkAffty(Pts,afftyPar)
%function [Pts,A,binNhbr,mdist] = mkAffty(Pts,afftyPar)
% I/P: Pts is an image array 
% O/P: Pts array has extra columns 
%
  
  %=========================
  %===== Parameters  =======
  %=========================
  dsThres = afftyPar.dsThres;
  dsSupp  = afftyPar.dsSupp;
  rho     = afftyPar.rho;
  sizeIm  = afftyPar.sizeIm;
  
  n = size(Pts);
  [xi, yi] = meshgrid(1:sizeIm(2), 1:sizeIm(1));
  Pts(:,2) = xi(:);
  Pts(:,3) = yi(:);
  
  %% Distance between consecutive points, and median distance
  X = Pts(:,1) * ones(1, n(1));
  X = X - X';
  Y = Pts(:,2) * ones(1, n(1));
  Y = Y - Y';
  Z = Pts(:,3) * ones(1, n(1));
  Z = Z - Z';
  
  dS  = max(abs(Y),abs(Z));
  d = (dS <= dsThres).* (X .^2);
  idx = find(d > 0);
  mdist = median(sqrt(d(idx)));
  sigma = rho*mdist; 
  binNhbr = (dS < dsThres) & (dS ~= 0);
  binSupp = (dS < dsSupp);
%  A = (dS <= dsThres).*exp( -d/(2 * sigma * sigma));
  minAffty = 0.01;
  A = (dS <= dsThres).*(exp( -d/(2 * sigma * sigma))+minAffty);
  A = sparse(A);

