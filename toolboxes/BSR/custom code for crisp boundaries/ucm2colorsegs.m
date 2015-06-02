% function [colorsegs] = ucm2colorsegs(ucm,im,k)
%
% Phillip Isola
% February 2014

function [colorsegs, labels] = ucm2colorsegs(ucm,im,k)
    
    %
    im = im2double(im);
        
    % get superpixels at scale k without boundaries:
    labels = bwlabel(ucm <= k);
    %labels = bwareaopen(labels,20);
    % make higher resolution
    smallRegionMerging_tweaked;
    res_scale = 1;
    labels = imresize(labels,res_scale*size(im(:,:,1)),'nearest');
    m = labels==0;
    m = bwmorph(m,'thin');
    labels = inPaint(labels,labels~=0);
    labels(m)=0;
    
    ii = unique(labels);
    %ii = ii(2:end); % remove zero label

%     colorsegs = ones([size(labels),size(im,3)]);
%     for i=1:length(ii)
%         m = labels==ii(i);
%         for c=1:size(im,3)
%             tmp = imresize(im(:,:,c),res_scale);
%             colorsegs_tmp = colorsegs(:,:,c);
%             colorsegs_tmp(m) = mean(tmp(m));
%             colorsegs(:,:,c) = colorsegs_tmp;
%         end
%     end  
    
    colorsegs = zeros([size(labels),size(im,3)]);
    for i = 1:length(ii)
        colorsegs = colorsegs + cat(3,S(ii(i)).color(1) * (labels == ii(i)),...
            S(ii(i)).color(2) * (labels == ii(i)), S(ii(i)).color(3) * (labels == ii(i)));
    end    
end