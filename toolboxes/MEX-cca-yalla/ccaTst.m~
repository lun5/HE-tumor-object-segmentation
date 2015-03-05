% CCA: connected component analysis

%==============================================
% connected components on the face image

FALSE = (0 == 1);
TRUE = ~FALSE;
%im = pgmRead('/u/jepson/tmp/chakraSmall.pgm');
fNameHdr = 'chakraSmall';
fName = [fNameHdr '.pgm'];
im = pgmRead(fName);

sizeIm = size(im);
useImage = TRUE;
Pts = ones(prod(sizeIm),2);
Pts(:,1) = im(:);


figure(502); clf;
showIm(reshape(Pts(:,1),sizeIm));

%%%  Affinity Matrix
rho = 1;
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

dsThres = 2;
dS  = max(abs(Y),abs(Z));
d = (dS < dsThres).* (X .^2);
idx = find(d > 0);
mdist = median(sqrt(d(idx)))
sigma = rho*mdist; % for face: 1
binNhbr = (dS < dsThres) & (dS ~= 0);
A = (dS < dsThres).*exp( -d/(2 * sigma * sigma));
imStats(A)

% A = affinity matrix, tau = threshold
tau = 0.6;
[lbl z] = cca(A,tau);


% each pixel is given a label stored in: lbl
ulbl = unique(lbl);
% total number of connected components
% the output variable z carries the same information
length(ulbl)

% sort the label information
for t = 1:length(ulbl)
  ulblLen(t) = length(find(lbl == ulbl(t)));
end
[sulblLen,iulblLen] = sort(-ulblLen);
sulblLen = -sulblLen;

% display connected components
kMin = 1; kMax = length(unique(lbl)); k = 1;
stop = FALSE; cnt = 1;
hFig = figure(502); clf;
while ~stop

  %ids = find(lbl==silbLen(k));
  ids = find(lbl==ulbl(iulblLen(k)));
  if (ids)
    figure(hFig);cla;
    %resizeImageFig(hFig, 0.5*size(im), 2); hold on;
    %set(get(hFig,'CurrentAxes'),'Ydir','reverse');
    showIm(reshape(im(:).*(lbl==iulblLen(k)),sizeIm));
    title(sprintf('Connected Component: %d',k));
  end

  [k stop cnt] = backForthButtonPress(hFig, k, kMin, kMax, cnt);
end



%==============================================
zz = find(lbl==1);

B = A(zz,zz);

Bd = sum(B,2);

% ceiling
thresh = 50;
tt =  find(Bd>=thresh);

% dumb loop
for i = 1:length(tt)
  sB = sum(B(:,tt(i)));
  nbrs = find(B(:,tt(i)));
  B(nbrs,tt(i)) = thresh/sB;
  B(tt(i),nbrs) = thresh/sB;    
end

%==============================================