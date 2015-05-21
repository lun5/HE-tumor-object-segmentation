function [n, r, converged, nVals, rVals, sumWs]=...
    robustIteration(n0, r0,  pEdgel, sigmaRho, nIts)
%% [n, r, converged, nVals, rVals ] = robustIteration(n0, r0, pEdgel,...
%%                                                   sigmaRho, nIts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Robust fitting iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  global FALSE;
  global TRUE;
  foundLine = FALSE;
  minWght0 = 1;
  n = n0(:); r = r0; sumW=0;
  nVals = n; rVals = r; nPrev = n0; sumWs = [];
  for it=1:nIts
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute errors and robust weights of active edgels 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute error residuals
    err = pEdgel*n+r; 
    
    % Compute robust weights
    W = 2*sigmaRho^2./(sigmaRho^2 + err.^2).^2;

   % Sum of all weights contributing to current segment.
    prevSumW = sumW;
    sumW = max(eps, sum(W));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check for sufficient weight and convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    converged = FALSE;
    if (it >= 2)
      if prevSumW > sumW - 1.0e-3;
         % Converged
         converged = TRUE; 
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve weighted eigenvalue problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b=sum([W,W].* pEdgel)'/sumW;
    d = pEdgel - repmat(b', size(pEdgel,1),1);
    C=(([W,W].*d )'*d)/sumW;
    
    
    %% Solve eigenvalue problem
    [V E] = eig(C);
    E = diag(E);
    minIndex = (E == min(E));
    if (sum(minIndex)>1)
      minIndex(2) = 0;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute fitted line parameters, n, c and t.
    %% Switch sign of normal, if necessary, to match previous iteration.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = V(:, minIndex);
    if ( nPrev' * n < 0.0)
      n = -n;
    end
    r = -n' * b;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Prepare for next iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nVals = [nVals n];
    rVals = [rVals r];
    sumWs = [sumWs; sumW];
    nPrev = n;
    prevSumW = sumW;
    if converged
      break;
    end

  end  % End of robust estimation iteration loop

