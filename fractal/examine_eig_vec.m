% Sample script for generating arbitrarily large
% fractal images and computing their affinity/laplacian 
% matrices

%[im,segIm] = rbfFracImage;

% 1 Gig: 2^30 this should be 
% total # of pixels = (image rows) x (image columns)
% => [rows cols] = sizeIm = [2^15 2^15];

% but on my 64-bit m/c i have to keep it small
% how much memory is there on this m/c anyways?
sizeIm = [128 128];
[im,segIm] = rbfFracImageNew(round(sum(1000*clock)), 2, 0.8, sizeIm);

figure(1); clf; showIm(im); figure(2); clf; showIm(double(segIm));
im = (im - min(im(:)))./(max(im(:)) - min(im(:)))*255;
im = uint8(im);

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

D = full(sum(A, 1)');      
sqrtD = D .^ 0.5; 
sqrtDinv = spdiags(sqrtD .^ -1,0,length(sqrtD),length(sqrtD));
L = sqrtDinv*A*sqrtDinv;

[U,S,V] = svds(L,30);
figure;
for i = 1:9
    tight_subplot(3,3,i); imagesc(reshape(U(:,i),size(im,1),size(im,2)));
    axis equal; 
    axis off;
    colormap('gray')
end

figure;
ha = tight_subplot(2,4,[.01 .0],[0 0],[0 0]);
% for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end 
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
for i = 1:8
    axes(ha(i));
    %subplot(3,3,i);
%     imagesc(reshape(EigVect(:,i), orig_sz(1), orig_sz(2)))
    imagesc(vect(:,:,i));
    axis equal; axis off; axis tight; colormap('gray');
end
tightfig;
%==============================================================
% in the dir: ~/ComputerVision/newMinCut/minCut/
clear all
close all

%im = imread('fracTest2.pgm');
[im,segIm] = rbfFracImageNew(round(sum(1000*clock)), 2, 0.8, [4 4]);
figure; imagesc(im); colormap('gray');
% rescale to be between 0 and 255
im = (im - min(im(:)))./(max(im(:)) - min(im(:)))*255;
im = uint8(im);
figure; imagesc(im); colormap('gray');
sizeIm = size(im);
Pts = ones(prod(sizeIm),2);
Pts(:,1) = im(:);

mDistMin = 2;  % mkAffinity does not care 
dsThres = 1; % distance threshold
afftyPar.dsThres = dsThres;
afftyPar.rho     = 1.5; % needed in mkAffinity
afftyPar.mdistMin     = mDistMin; % mkAffinity does not care
afftyPar.dsSupp = 0;
afftyPar.sizeIm = sizeIm;
afftyPar

[Pts,A,binNhbr,binSupp,mdist] = mkAffty(Pts,afftyPar);
%A = brightAfftyNew(im,dsThres,max(mDistMin,10));% Fast!
B = brightAfftyNew(im,dsThres,mdist*afftyPar.rho);% Fast!
%figure; spy(A);

D = full(sum(A, 1)');      
sqrtD = D .^ 0.5; 
sqrtDinv = spdiags(sqrtD .^ -1,0,length(sqrtD),length(sqrtD));
L = sparse(sqrtDinv*A*sqrtDinv);

[U,S,V] = svds(L,20);
figure;
for i = 1:9
    subplot(3,3,i); imagesc(reshape(U(:,i),sizeIm(1),sizeIm(2)),[min(U(:)) max(U(:))]);
    axis off;
    colormap('gray')
end

figure;
for i = 10:18
    subplot(3,3,i-9); imagesc(reshape(U(:,i),sizeIm(1),sizeIm(2)));
    axis off;
    colormap('gray')
end

%==============================================================
% Luong's code

% 
% %% Normalized graph laplacian 
% [wx, wy] = size(W);
% x = 1 : wx;
% S = full(sum(W, 1));
% D = sparse(x, x, S, wx, wy);
% %clear S x;
% 
% opts.issym=1;opts.isreal = 1;opts.disp=0;
% nvec = min(30,size(D,1));
% % generalized eigen problem (D-W)u = lambda Du
% % equivalent to inv(D) (D-W) u = lambda u
% [EigVect, EVal] = eigs(D - W, D, nvec, 'sm',opts);
% %clear D W opts;
% 
% EigVal = diag(EVal);
% %clear EVal;
% % arrange in ascending order
% EigVal(1:end) = EigVal(end:-1:1);
% EigVect(:, 1:end) = EigVect(:, end:-1:1);
% 
% %%
% txo=size(im,1); tyo=size(im,2);
% vect = zeros(txo, tyo, nvec);
% for v = 2 : nvec,
%     vect(:, :, v) = reshape(EigVect(:, v), [tyo txo])';
% end
% %montage2(vect);
% %clear EigVect;
% 
% 
% figure;
% for i = 1:9
%     subplot(3,3,i); imagesc(vect(:,:,i));
%     axis equal;axis([0 tyo 0 txo])
%     set(gca,'xtick',[]);set(gca,'ytick',[]);
%     colormap('gray')
% end


