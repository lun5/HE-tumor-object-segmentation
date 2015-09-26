function [classes,X,lambda,Xr,W,C,timing] = ncut_multiscale(image,nsegs,options);
% [classes,X] = ncut_multiscale(image,10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                  %
%   Multiscale Normalized Cuts Segmentation Code   %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% inputs: 
% image: image to segment (size pxq or pxqx3)
% nsegs: number of segments requested
% outputs:
% classes: image regions (size pxq)
% X: eigenvectors (size pxqxnsegs)
% lamda: eigenvalues
% Xr: rotated eigenvectors (computed during discretisation)
% W: multiscale affinity matrix
% C: multiscale constraint matrix
% timing: timing information
% 
% source code available at http://www.seas.upenn.edu/~timothee
% Authors: Timothee Cour, Florence Benezit, Jianbo Shi
% Related publication:
% Timothee Cour, Florence Benezit, Jianbo Shi. Spectral Segmentation with
% Multiscale Graph Decomposition. IEEE International Conference on Computer
% Vision and Pattern Recognition (CVPR), 2005.
% 
% Please cite the paper and source code if you are using it in your work.


image=im2double(image);
[p,q,r] = size(image);

if nargin<3
    options=[];
end

%% compute multiscale affinity matrix W and multiscale constraint matrix C
t= cputime;
[layers,C]=compute_layers_C_multiscale(p,q);
dataW = computeParametersW(image);
W=computeMultiscaleW(image,layers,dataW,options);
disp(['time compute W,C: ',num2str(cputime-t)]);

%% compute constrained normalized cuts eigenvectors
disp('starting multiscale ncut...');
t = cputime;
if ~isempty(C)
    [X,lambda,timing] = computeNcutConstraint_projection(W,C,nsegs);
else
    [X,lambda,timing] = computeKFirstEigenvectors(W,nsegs);
end

timing.total=cputime-t;

disp(['time ncut: ',num2str(timing.total)]);
% disp(['time W*X: ',num2str(timing.A_times_X)]);
% disp(['time dsaupd: ',num2str(timing.dsaupd)]);
% disp(['nbA_times_X: ',num2str(timing.nbA_times_X)]);


%% compute discretisation
[p,q,r]=size(image);
indPixels = (1:p*q)';
X = reshape(X(indPixels,:),p,q,nsegs);
t =cputime;
[classes,Xr] =discretisation(X);
disp(['time discretize: ',num2str(cputime-t)]);

