function [HE, H, E] = extractStainPixels(I, Io, beta, alpha)
% Extract example of purple and pink pixels from 
% Input:
% I         - RGB input image;
% Io        - (optional) transmitted light intensity (default: 240);
% beta      - (optional) OD threshold for transparent pixels (default: 0.15);
% alpha     - (optional) tolerance for the pseudo-min and pseudo-max (default: 1);
% HERef     - (optional) reference H&E OD matrix (default value is defined);
%
% Output:
% HE        - stain vectors in OD; H first column, E second column
% H         - hematoxylin pixel (RGB);
% E         - eosin pixel (RGB);
%
% References:
% A method for normalizing histology slides for quantitative analysis. M.
% Macenko et al., ISBI 2009
%
% transmitted light intensity
% Luong Nguyen 9/14/15
% Inspired by code written by Mitko Veta
if ~exist('Io', 'var') || isempty(Io)
    Io = 240;
end

% OD threshold for transparent pixels
if ~exist('beta', 'var') || isempty(beta)
    beta = 0.15;
end

% tolerance for the pseudo-min and pseudo-max
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = [1 2];
end

I = double(I);
I = reshape(I, [], 3);

% calculate optical density
OD = -log((I+1)/Io);

% remove transparent pixels
ODhat = OD(~any(OD < beta, 2), :);
Ihat = I(~any(OD <beta,2),:);
%% calculate eigenvectors, why not SVD????
[V, ~] = eig(cov(ODhat));

% project on the plane spanned by the eigenvectors corresponding to the two
% largest eigenvalues
That = ODhat*V(:,2:3);

% find the min and max vectors and project back to OD space
phi = atan2(That(:,2), That(:,1));

minPhi = prctile(phi, alpha); 
maxPhi = prctile(phi, 100-alpha);

index_minPhi = phi >= minPhi(1) & phi <= minPhi(2);
index_maxPhi = phi >= maxPhi(2) & phi <= maxPhi(1);

vMin = V(:,2:3)*[cos(minPhi(1)); sin(minPhi(1))];
vMax = V(:,2:3)*[cos(maxPhi(1)); sin(maxPhi(1))];

% a heuristic to make the vector corresponding to hematoxylin first and the
% one corresponding to eosin second
if vMin(1) > vMax(1)
    HE = [vMin vMax];
    H = Ihat(index_minPhi,:); E = Ihat(index_maxPhi,:);
else
    HE = [vMax vMin];
    H = Ihat(index_maxPhi,:); E = Ihat(index_minPhi,:);
end

end
