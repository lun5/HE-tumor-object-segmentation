% Sample script for generating arbitrarily large
% fractal images and computing their affinity/laplacian 
% matrices

%[im,segIm] = rbfFracImage;

% 1 Gig: 2^30 this should be 
% total # of pixels = (image rows) x (image columns)
% => [rows cols] = sizeIm = [2^15 2^15];

% but on my 64-bit m/c i have to keep it small
% how much memory is there on this m/c anyways?
sizeIm = [64 64];
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

W = brightAfftyNew(im,2,max(mDistMin,10));% Fast!

%% Normalized graph laplacian 
[wx, wy] = size(W);
x = 1 : wx;
S = full(sum(W, 1));
D = sparse(x, x, S, wx, wy);
%clear S x;

opts.issym=1;opts.isreal = 1;opts.disp=0;
nvec = min(30,size(D,1));
% generalized eigen problem (D-W)u = lambda Du
% equivalent to inv(D) (D-W) u = lambda u
[EigVect, EVal] = eigs(D - W, D, nvec, 'sm',opts);
%clear D W opts;

EigVal = diag(EVal);
%clear EVal;
% arrange in ascending order
EigVal(1:end) = EigVal(end:-1:1);
EigVect(:, 1:end) = EigVect(:, end:-1:1);

%%
txo=size(im,1); tyo=size(im,2);
vect = zeros(txo, tyo, nvec);
for v = 2 : nvec,
    vect(:, :, v) = reshape(EigVect(:, v), [tyo txo])';
end
%montage2(vect);
%clear EigVect;

figure;
for i = 1:9
    subplot(3,3,i); imagesc(vect(:,:,i));
    axis equal;axis([0 tyo 0 txo])
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    colormap('gray')
end


