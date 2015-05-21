function plotLike(modelSize, fitResults, crossValResults, h)
% plotLike(modelSize, fitResults, crossValResults, h)
%
% PLOTLIKE: Plot log likelihood results for fitting 1D mixture models
% Example Useage: 
%    See mixScript.m
% Input:
%   modelSize:  range of nInliers eg. 0:10
%   fitResults(:,1) log likelihood fit results for range of inliers 
%   fitResults(:,2:3)  error bars on log like fit results.
%   crossValResults(:,1) data log likelihood on test set. (=[] if no cross
%                        validation was used.)
%   h: figure handle
 if nargin<4
   h = figure(2);
 end

figure(h);
 clf;
 hAxes = gca;
 set(hAxes, 'ColorOrder', [0 1 0], ...
    'LineStyleOrder', '-');
 set(hAxes,'FontSize', 16);
 
 if size(fitResults,2) > 1
   errorBars = 1;
   eBarX = zeros(size(fitResults,1),2);
   eBarY = eBarX;
   eBarX(:,1) = modelSize(:);
   eBarX(:,2) = modelSize(:);
   eBarY(:,1) = fitResults(:,2) + fitResults(:,3);
   eBarY(:,2) = fitResults(:,2) - fitResults(:,3);
 else
   errorBars = 0;
 end

 nCrossVal = size(crossValResults,2);
 if nCrossVal > 0
   [mL ML] = range2([fitResults(:,1) crossValResults(:,1:nCrossVal)]);
 else
   [mL ML] = range2(fitResults(:,1));
 end

 mL = floor(mL/50) * 50; 
 ML = ceil(ML/50) * 50;
 axis([0.0 10.5 mL ML]);
 hold on;
 if nCrossVal == 0
   hLegend = zeros(1+errorBars,1);
 else
   hLegend = zeros(2+errorBars,1);
 end
 kLegend = 1;
 hLegend(1) = plot(modelSize(:), fitResults(:,1), 'rx-', 'LineWidth',2);

 if errorBars
   plot(modelSize(:), fitResults(:,2), 'go');
   kLegend = kLegend+1;
   hLegend(kLegend) = plot(modelSize(:), fitResults(:,2), 'g');
   line(eBarX', eBarY');
 end
  
 xlabel('Number of Inlier Components'); ylabel('Log(Likelihood)');

 if nCrossVal > 0
   kLegend = kLegend+1;
   for i=1:nCrossVal
     hLegend(kLegend) = plot(modelSize(:), crossValResults(:,i), 'bo-', ...
                             'LineWidth',2);
   end
 end

 if kLegend == 3
   legend(hLegend,'Fit Model','Exp. Cross Validation',...
	  'Data Cross Validation',4);
   title('Fit Likelihood and Cross Validation');
 elseif kLegend == 2
   legend(hLegend,'Fit Model', 'Data Cross Validation',4);
   title('Fit Likelihood and Cross Validation');
 else
   legend(hLegend,'Fit Model',4);
   title('Fit Likelihood');
 end

 return
