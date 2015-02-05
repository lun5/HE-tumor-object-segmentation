

% $Author: ChrisMcCormick $    $Date: 2014/05/19 22:00:00 $    $Revision: 1.3 $
% Mixture of von Mises distribution with an uniform noise component
addpath(genpath(pwd));
addpath(genpath(fullfile(pwd,'toolboxes','CircStat2012a')));
%%======================================================
%% STEP 1a: Generate data from two 1D distributions.
% clear vars; close all;
% mu1 = 0;      % Mean
% kappa1 = 70;    % kappa
% m1 = 300;      % Number of points
% 
% mu2 = pi;
% kappa2 = 10;
% m2 = 500;
% 
% mu3 = 5*pi/4;
% kappa3 = 0;
% m3 = 100;
% 
% % Generate the data.
% X1 = circ_vmrnd(mu1,kappa1, m1);
% X2 = circ_vmrnd(mu2,kappa2, m2);
% X3 = circ_vmrnd(mu3,kappa3, m3);
% 
% % have to convert this into cartesian coordiates
% X1_cart = [cos(X1) sin(X1)];
% X2_cart = [cos(X2) sin(X2)];
% X3_cart = [cos(X3) sin(X3)];
% 
% X = [X1; X2; X3];
% X_cart = [X1_cart; X2_cart;X3_cart]; 
% d = size(X_cart,2); % dimension
% %%=====================================================
% %% STEP 1b: Plot the data points and their pdfs.
% 
% x = -pi:0.1:pi;
% y1 = circ_vmpdf(x, mu1, kappa1);
% y2 = circ_vmpdf(x, mu2, kappa2);
% y3 = circ_vmpdf(x, mu3, kappa3);
% 
% figure;
% plot(x, y1, 'b-');
% hold on;
% plot(x, y2, 'r-');
% plot(x, y3, 'g-');
% plot(X1, zeros(size(X1)), 'bx', 'markersize', 10);
% plot(X2, zeros(size(X2)), 'rx', 'markersize', 10);
% plot(X3, zeros(size(X3)), 'g+', 'markersize', 10);
% 
% xlim([-pi pi]);
% set(gcf,'color','white') % White background for the figure.
% hold off

%% test with the image
datadir = 'test_images'; 
opts_input = setEnvironment_inputs;
imname = 'gland3.tif';
I = getImage(datadir, imname, opts_input);

thetas = load('theta_test.mat');
thetas = thetas.f_maps{1};
X = thetas(:);
X_cart = [cos(X) sin(X)];

k = 3;
%% Call the function
opts.noise = 1;
[ mu_hat_polar,mu_hat_cart, kappa_hat,posterior_probs, prior_probs] = ...
     moVM(X_cart,k,opts);
[idxbest, centroids_cart] = spkmeans(X_cart,k,opts);
% membership
[~, indx_membership] = max(posterior_probs,[],2); % 4 is the uniform noise

id1 = reshape(indx_membership, size(thetas));
id1(id1 ~=1) = 0;
id1_im = I.*uint8(repmat(id1,1,1,3));
figure; imshow(id1_im);
set(gcf,'color','white') % White background for the figure.

id2 = reshape(indx_membership, size(thetas));
id2(id2 ~=2) = 0;
id2(id2 == 2) = 1;
id2_im = I.*uint8(repmat(id2,1,1,3));
figure; imshow(id2_im);
set(gcf,'color','white') % White background for the figure.

id3 = reshape(indx_membership, size(thetas));
id3(id3 ~=3) = 0;
id3(id3 == 3) = 1;
id3_im = I.*uint8(repmat(id3,1,1,3));
figure; imshow(id3_im);
set(gcf,'color','white') % White background for the figure.

id4 = reshape(indx_membership, size(thetas));
id4(id4 ~=4) = 0;
id4(id4 == 4) = 1;
id4_im = I.*uint8(repmat(id4,1,1,3));
figure; imshow(id4_im);
set(gcf,'color','white') % White background for the figure.

%%=====================================================
%% Plot the data points and their estimated pdfs.

x = -pi:0.1:pi;
c = ['r','g','b'];

%figure;
% for i=1:k
%     yk = circ_vmpdf(x, mu_hat_polar(i), kappa_hat(i));
%     plot(x, yk,'Color',c(i),'LineStyle','-'); hold on;
%     plot(x, yk,'Color','k','LineStyle','-'); 
% end
% hold off

figure;
for i=1:k
    yk = circ_vmpdf(x, mu_hat_polar(i), kappa_hat(i));
    plot(x, yk,'Color',c(i),'LineStyle','-'); hold on;
    %plot(x, yk,'Color','k','LineStyle','-'); 
end
histogram(X,'Normalization','pdf','FaceColor',[0.8 0.8 0.8]);

hold off
xlim([-pi pi]);
set(gcf,'color','white') % White background for the figure.
%hold off
% Plot over the existing figure, using black lines for the estimated pdfs.
% figure;
% plot(x, y1, 'k-'); hold on
% plot(x, y2, 'k-');
% plot(x,y3, 'k-'); 
% hold off;

%[~,cluster_id] = max(posterior_probs,[],2); 
%groundtruth_indx = [ones(m1,1)*2; ones(m2,1)];

