%% 
addpath(genpath(pwd));
% I = imread('test_images/253027.jpg'); % zebra
% I = imread('test_images/syntheticImage.png'); % synthetic image
% I = imread('test_images/tp09-96_20480_10240_2048_2048.tif'); % H&E
%datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
%I = imread(fullfile(datadir,'tp10-867-1_47104_22528_2048_2048.tif'));
%I = imread(fullfile(datadir,'tp10-867-1_26624_24576_2048_2048.tif'));
%I = imread(fullfile(datadir,'tp10-867-1_34816_18432_2048_2048.tif'));
%I = imread(fullfile(datadir,'tp10-611_22528_16384_2048_2048.tif'));
% imshow(I);
% rect = getrect;
% I = imcrop(I,rect);
opts = setEnvironment('speedy');

datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
if ~ exist(datadir,'dir')
    datadir = '/Users/lun5/Research/color_deconvolution/TissueImages/';
end
I = imread(fullfile(datadir,'tp10-867-1_4096_20480_2048_2048.tif'));
rect = [440         746        1178         489];
I = imcrop(I,round(rect));

I = im2uint8(I);
if (size(I,3)==1)
    I = repmat(I,[1 1 3]);
end

num_scales = opts.num_scales;
scale_offset = opts.scale_offset;
f_maps = getFeatures(double(I)/255,num_scales+scale_offset,opts.features.which_features,opts);
%getFeatures_theta;
Nsamples = 10000;%.kde.Nkernels;
%opt.sig = 20;
F = sampleF(f_maps,Nsamples,opts);
Fsym = [F; [F(:,2) F(:,1)]]; % symmetric F(A,B) = F(B,A). 
p = kde(Fsym',0.05,[],'e');
    
 figure; ndhist(F(:,1),F(:,2),'axis',[0 1 0 1],'filter','bins',1,'columns');
 xlabel('Luminance A'); ylabel('Luminance B');
 %figure; ndhist(F(:,1),F(:,2),'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
 %xlabel('\phi'); ylabel('\psi'); 
 axis square;set(gca,'FontSize',16);set(gcf,'color','white') 
 
 % PMI
 tol = opts.kde.kdtree_tol;
 pJoint = evaluate(p,F',tol);
 prct = 10;
 reg = prctile(nonzeros(pJoint),prct); %opts.p_reg;
 pJoint = pJoint + reg;
 N = floor(size(F,2)/2); assert((round(N)-N)==0);
 p2_1 = marginal(p,1:N);
 p2_2 = marginal(p,N+1:(2*N));
 p2 = joinTrees(p2_1,p2_2,0.5);
 % pMarg_x = evaluate_batches(p2,X(:)',tol);
 % pMarg_y = evaluate_batches(p2,Y(:)',tol);
 pMarg_x = evaluate(p2,F(:,1)',tol);
 pMarg_y = evaluate(p2,F(:,2)',tol);
 pProd = pMarg_x.*pMarg_y +reg;

 %% calculate pmi
 pmi = (pJoint.^(opts.joint_exponent))./pProd;
 log_pmi = log(pmi);
 
 %% absolute distances
 minAffty = 0.01;
 Fdist = sqrt((F(:,1)-F(:,2)).^2);
 mDist = median(Fdist);
 Fdist =(exp(-(Fdist.^2)/(2*mDist^2))+minAffty);
 %
 figure; plot(Fdist, log_pmi,'ro','MarkerSize',3);
 xlabel('Difference Measure'); ylabel('logPMI');
 xlim([0 1]);set(gca,'FontSize',16);set(gcf,'color','white') 
 %
 figure; plot(Fdist, pmi,'k+','MarkerSize',3);
 xlabel('Difference Measure'); ylabel('PMI');
 xlim([0 1]);set(gca,'FontSize',16);set(gcf,'color','white') 
 %
%  figure; plot(blockSizes,times,'k-','LineWidth',3);
%  xlabel('Block Sizes (pixel)'); ylabel('Time (s)');
%  set(gca,'XTick',blockSizes); 
%  xlim([128 2048]);
%  
 
%% create the contour plots for P_AB
% look up help kde/hist
tol = opts.kde.kdtree_tol;
x = 0:0.01:1;y = 0:0.01:1;
%thetaRange = pi;
%x = -thetaRange:0.01:thetaRange;y = -thetaRange:0.01:thetaRange;

[X,Y] = meshgrid(x,y);
Fim = [X(:),Y(:)];
pd = evaluate(p,Fim',tol);
pd_mesh =  reshape(pd, size(X));
%[pd,x,y] = hist(p,500,[1,2]);
% 
figure;mesh(x,y,pd_mesh); axis square; colorbar;
xlabel('Luminance A'); ylabel('Luminance B');
set(gca,'XTick',0:0.1:1)
set(gca,'YTick',0:0.1:1)


figure;contourf(x,y,pd_mesh,20); axis square; colorbar;
xlabel('Luminance A'); ylabel('Luminance B');
set(gca,'XTick',0:0.1:1);set(gca,'YTick',0:0.1:1)
%set(gca,'FontSize',30); set(gca,'XTick',-pi:pi/2:pi,'YTick', -pi:pi/2:pi);

%% Interactive selection of red, green, blue circles on the zebra
figure; imshow(I); hold on;

[rowSub,colSub] = ginput;
rowSub = round(rowSub); colSub = round(colSub); 
% c_vecs = {'r','r','g','g','w','w'};
c_vecs = {'r','g','b'};
% shape inserter for 6 combinations: 
% pink-pink: red circle, purple purple: white circle, white white: green
% circle, pink purple: white square, pink white: red square, purle-white:
% green square

for i =1:floor(length(colSub)/2)
    coord1 = [rowSub(2*(i-1)+1),colSub(2*(i-1)+1)];
    coord2 = [rowSub(2*i),colSub(2*i)];
    r_shape = 2; % radius of the shape drawn
    plot([coord1(1),coord2(1)],[coord1(2),coord2(2)],'o','MarkerSize',4,...
        'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
%     if mod(i,2) == 1
%         viscircles((coord1 + coord2)/2,r_shape, 'EdgeColor',c_vecs{i},'DrawBackgroundCircle',false);
%     else       
        rectangle('Position',[(coord1(1) + coord2(1))/2 - r_shape ...
            (coord1(2) + coord2(2))/2 - r_shape 2*r_shape 2*r_shape],...,
            'LineWidth', 2, 'EdgeColor',c_vecs{i})
%     end
end

hold off;

%% 
%linearInd = sub2ind(size(f_maps{1}), rowSub, colSub);
% if later not work, add {1} behind f_maps
f_maps_cur = f_maps{1};
linearInd = sub2ind(size(f_maps_cur), colSub, rowSub);
feats = f_maps_cur(linearInd);

%% joint probabilities
prct = 5;
reg = prctile(nonzeros(pd),prct); %opts.p_reg;
pJoint = reg + pd;
% in a mesh
pJoint_mesh = reg + pd_mesh;%reshape(pJoint, size(X));
%pJoint_mesh = log(pJoint_mesh);
%pJoint_mesh = reg + pd; 

figure;mesh(x,y, pJoint_mesh); axis square; colorbar;
xlabel('Luminance A'); ylabel('Luminance B');
set(gca,'XTick',0:0.2:1);set(gca,'YTick',0:0.2:1)
set(gca,'FontSize',16);


figure;contourf(x,y,pJoint_mesh,20); axis square; colorbar;
xlabel('Luminance A'); ylabel('Luminance B');
%xlabel('Theta A'); ylabel('Theta B');
set(gca,'XTick',0:0.2:1);set(gca,'YTick',0:0.2:1)
set(gca,'FontSize',16);
hold on;
for i =1:floor(length(feats)/2)
    coord = [feats(2*(i-1)+1) feats(2*i)]; % coordinate of points picked interactively
    r_shape = 0.03;
    plot(coord(1),coord(2),'o','MarkerSize',4,...
        'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
%     if mod(i,2) == 1
%         viscircles(coord,r_shape, 'EdgeColor',c_vecs{i},'DrawBackgroundCircle',false);
%     else 
        rectangle('Position', [ coord(1) - r_shape coord(2) - r_shape ...
            2*r_shape 2*r_shape], 'LineWidth', 2, 'EdgeColor',c_vecs{i});
%     end
end
hold off;
%% useful functions to look up
% help kde/evaluate
% help kde/getPoints
%% evaluate p(A)p(B)
N = floor(size(Fim,2)/2); assert((round(N)-N)==0);
p2_1 = marginal(p,1:N);
p2_2 = marginal(p,N+1:(2*N));
p2 = joinTrees(p2_1,p2_2,0.5);
% pMarg_x = evaluate_batches(p2,X(:)',tol);
% pMarg_y = evaluate_batches(p2,Y(:)',tol);
pMarg_x = evaluate(p2,X(:)',tol);
pMarg_y = evaluate(p2,Y(:)',tol);
pProd = pMarg_x.*pMarg_y +reg;

%% calculate pmi
%log_pmi = log((pJoint.^(opts.joint_exponent))./pProd);
pmi = (pJoint.^(opts.joint_exponent))./pProd;
logpmi = log(pmi);
normalizedPMI = logpmi./(-log(pJoint)); 
Z_pmi = reshape(logpmi,size(X));
%Z_pmi = reshape(log(pmi),size(X));
%Z_pmi = griddata(x,y,log_pmi,X,Y,'cubic');
% 
figure;mesh(x,y,Z_pmi);axis square; colorbar;
xlabel('Luminance A'); ylabel('Luminance B');
set(gca,'XTick',0:0.2:1);set(gca,'YTick',0:0.2:1)
set(gca,'FontSize',16);

% 
% 
figure;[C_pmi,h_pmi]=contourf(x,y,Z_pmi,20);
axis square; colorbar; 
xlabel('Luminance A'); ylabel('Luminance B');
% xlabel('Theta A'); ylabel('Theta B');
set(gca,'XTick',0:0.2:1);set(gca,'YTick',0:0.2:1)
set(gca,'FontSize',16);

hold on;
for i =1:floor(length(feats)/2)
    coord = [feats(2*(i-1)+1) feats(2*i)]; % coordinate of points picked interactively
    r_shape = 0.03;
    plot(coord(1),coord(2),'o','MarkerSize',4,...
        'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
%     if mod(i,2) == 1
%         viscircles(coord,r_shape, 'EdgeColor',c_vecs{i},'DrawBackgroundCircle',false);
%     else 
        rectangle('Position', [ coord(1) - r_shape coord(2) - r_shape ...
            2*r_shape 2*r_shape], 'LineWidth', 2, 'EdgeColor',c_vecs{i});
%     end
end
hold off; 

%close all;

