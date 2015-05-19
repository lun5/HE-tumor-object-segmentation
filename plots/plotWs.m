%coords = [114,45;71,35;168,87;169,44]; %leopard
% imshow(I);
% [rowSub,colSub] = ginput;
% coords = round([rowSub,colSub]);
%% Find the coordinates of the centers, every 2rad + 1 away 
rad = opts.localPairs.rad;
dist = round(rad/4); %rad
rowSub = rad+dist:(2*rad+dist):size(I,1);
colSub = rad+dist:(2*rad+dist):size(I,2);
[XX, YY] = meshgrid(rowSub,colSub);
coords = round([YY(:),XX(:)]);
    
Ws_current = Ws{opts.num_scales};

% figure; imshow(I./(2*255)); hold on;
% plot(coords(:,1), coords(:,2),'wx','MarkerSize',10);
% 
% hold off;
ind_coords = sub2ind([size(I,2),size(I,1)],coords(:,1),coords(:,2));
Ws_spots = Ws_current(ind_coords,:);
aff_im = zeros(size(I(:,:,1)));
for i = 1:3:size(coords,1)
    im = reshape(Ws_spots(i,:),size(I,2), size(I,1));
    aff_im = aff_im + im';
end
% movie(M,5,3);
% movie2avi(M,fullfile(pwd,'results','opp.avi'), 'compression', 'None','fps',3);
% movie2avi(M,fullfile(pwd,'results','pink_pink.avi'), 'compression', 'None','fps',3);

figure; imshow(I./255); hold on; 
%cspy(min(aff_im,15),'markersize',5, 'colormap', 'jet', 'levels', 7); 
cspy(min(aff_im,prctile(aff_im(:),90)),'markersize',5, 'colormap', 'jet', 'levels', 7);
%cspy(aff_im,'markersize',5, 'colormap', 'jet', 'levels', 7);
colorbar; 
hold off;

% figure;
% h = histogram(aff_im(aff_im > 0),'DisplayStyle','bar' );
% h.FaceColor = [0.8 .8 .8];h.BinWidth = max(nonzeros(Ws_spots))/50;
% h.Normalization = 'probability'; ax = axis;
% axis([min(nonzeros(Ws_spots)) max(nonzeros(Ws_spots)) 0 1]);

Ws{opts.num_scales}= min(Ws{opts.num_scales},prctile(aff_im(:),90));