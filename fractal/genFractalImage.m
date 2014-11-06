% Sample script for generating arbitrarily large
% fractal images and computing their affinity/laplacian 
% matrices

%[im,segIm] = rbfFracImage;

% 1 Gig: 2^30 this should be 
% total # of pixels = (image rows) x (image columns)
% => [rows cols] = sizeIm = [2^15 2^15];

% but on my 64-bit m/c i have to keep it small
% how much memory is there on this m/c anyways?
sizeIm = [2^8 2^8];
[im,segIm] = rbfFracImageNew(round(sum(1000*clock)), 2, 0.8, sizeIm);

figure(1); clf; showIm(im); figure(2); clf; showIm(double(segIm));

% Build affinity and the laplacian matrices
mDistMin=2;
%spectralStep2;

FALSE = (0 == 1);
TRUE = ~FALSE;
MAX_UINT = 2^31;% on a 32-bit machine

sizeIm = size(im);
mask = FALSE;
fprintf('Building affinity matrix...');
%=====================================================
%========== affinity matrix            ===============
%=====================================================
afftyPar.dsThres = 1.1;
afftyPar.rho     = 1.5; 
afftyPar.mdistMin     = mDistMin; 
afftyPar

A = brightAfftyNew(im,2,max(mDistMin,10));% Fast!

