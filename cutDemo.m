addpath('data');

clear all
close all

FALSE = (0 == 1);
TRUE = ~FALSE;

%=====================================================
%========== input image selection      ===============
%=====================================================
[im,fNameHdr] = imForCutdemo;
sizeIm = size(im);
Pts = ones(prod(sizeIm),2);
Pts(:,1) = im(:);

clc;
fprintf('Building affinity matrix...');
%=====================================================
%========== affinity matrix            ===============
%=====================================================
afftyPar.sizeIm  = sizeIm;
afftyPar.dsThres = 1.1; %kanisza
afftyPar.rho     = 1.5; %kanisza
if (strcmp(fNameHdr,'chakraSmall') ==1 )
  afftyPar.rho     = 1.0; 
end
afftyPar

% Pts array is updated.
[Pts,A,binNhbr,mdist] = mkAffty(Pts,afftyPar);

fprintf('...Done \n');
fprintf('Hit enter to continue\n');
pause;
clc;
%=====================================================
%========== define cut parameters      ===============
%=====================================================
% only parameter that needs to be set for iterative 
% cutting
cutPar.beta0 = 40; % kanisza (20x24): 20
                   % kanisza (40x48): 80

cutPar.half0          = cutPar.beta0/4;  % /10  or half0 = 10;
cutPar.dLogHalfMin    = -0.2;            % dont change this
cutPar.dLogHalfMaxCut = -0.02;           % dont change this

cutPar.maxIts   = 1000;        % dont change this
cutPar.maxBasis = size(Pts,1); % fixed
cutPar.its      = 0; % set by itCut3
cutPar.dim      = 0; % set by itCut3

cutPar.sizeIm      = sizeIm;
cutPar.svdDone     = FALSE;
cutPar.useImage    = TRUE;
cutPar.displayON   = 1; % set this to 0 for no display
cutPar.displayStep = 1; % which iterations to show

clc

fprintf('EigenCut...\n');
%=======================================================
%========== itCut as a routine          ================
%=======================================================
cutPar
binCut0 = zeros(prod(sizeIm), prod(sizeIm));
[binCut0,cutPar,discComp,U,S,V] = itCut3(Pts,A,binCut0,binNhbr,cutPar);
displayComp(Pts,discComp,cutPar.useImage,cutPar.sizeIm);
fprintf('...EigenCut done\n');
clc

%=======================================================
%========== components: paint-by-numbers/colors  =======
%=======================================================
fprintf('Paint Components...\n');
ok = 1;
while(ok==1)
  displayComp(Pts,discComp,cutPar.useImage,cutPar.sizeIm);
  ok = input(sprintf('Repeat 1/0: '));
end
fprintf('...painting components done\n');
clc;

%=======================================================
%========== Markov kernel diffusion       ==============
%=======================================================
fprintf('Markov Kernel Diffusion...\n');
Acut = A .* ~binCut0; 
D = sum(Acut, 1)';              % Normalize column sum to one.
sqrtD = D .^ 0.5;
ok = 1;
mnItr = cutPar.dim;
mxItr = prod(cutPar.sizeIm);
while (ok==1)
  itr = input(sprintf('diffusion kernel power between(0-inf): '));
  %itr = max(itr,mnItr);
  %itr = min(itr,mxItr);
  [blurIm,K] = aDiffuse(im,U,S,V,sqrtD,itr,cutPar.sizeIm);
  ok = input(sprintf('Repeat 1/0: '));
end
fprintf('...kernel diffusion done\n');

%=======================================================
%========== worth saving results?       ================
%=======================================================
% save image
% save iCut results


