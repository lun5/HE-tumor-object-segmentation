function entropy = entMix1D(mixDist)
% entropy = entMix1D(mixDist);
% Estimates integral [- mixDist(x) log_e(mixDist(x))]dx

  entropy = crossEntMix1D(mixDist, mixDist);
  return;
