%% function [f_maps] = getFeatures(im_rgb,scale,which_feature,opts)
% 
% INPUTS
%  im_rgb        - NxMx3 query image; C can be 1 (grayscale) or 3 (rgb color)
%  scale         - how many times should the image be downsampled?
%  which_feature - which feature types to compute? can contain multiple entries
%  opts          - parameter settings (see setEnvironment)
%
% OUTPUTS
%  f_maps        - NxMxF array of F feature maps input im_rgb
% 
% -------------------------------------------------------------------------
% HE segmentation toolbox
% Luong Nguyen, 2014 [lun5@pitt.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

function [f_maps] = getFeatures(im_rgb,scale,which_features,opts)
    
    im = [];
    if isempty(opts);
        opts = setEnvironment_affinity;
    end
    %%
    %if strcmp(which_feature,'hue opp') || (strcmp(which_feature,'saturation opp'))  || (strcmp(which_feature,'brightness opp'))
    if ismember('hue opp',which_features) || (ismember('saturation opp',which_features))  || (ismember('brightness opp',which_features))
        %tmp = load(fullfile(pwd,'DanTrainingData','rotation_matrix_after_centering.mat'));
        %opts.features.rotation_matrix = tmp.rotation_mat;
        %opts.features.decorrelate = 0;
        r = im_rgb(:,:,1); 
        X = reshape(im_rgb,[size(im_rgb,1)*size(im_rgb,2),size(im_rgb,3)]);
        rotated_coordinates = opts.features.rotation_matrix*X'; %double([r(:)'; g(:)'; b(:)']);
        white_index = rotated_coordinates(1,:) > sum(opts.features.rotation_matrix(1,:)) - 1e-3;
    end
    
    for feature_iter = 1: length(which_features)
        
    if (strcmp(which_features{feature_iter},'luminance'))
        im = mat2gray(mean(im_rgb,3));
        if opts.affinity.plot; figure; imshow(im); end %image(im); axis off; axis equal; 
    end
    if (strcmp(which_features{feature_iter},'color'))
        %% color features
        if (size(im_rgb,3)==3)
            colorTransform = makecform('srgb2lab');
            cf = mat2gray(applycform(im_rgb,colorTransform));
            im = cat(3,im,cf(:,:,1:3));
        elseif (size(im_rgb,3)==1) % grayscale image
            im = im_rgb;
        else
            error('unhandled image format');
        end
        % NEED TO RESCALE THIS TO [0 255] 
    end
    
    if (strcmp(which_features{feature_iter},'var'))    
        %% variance features
        f = pcaIm(im_rgb);      
        Nhood_rad = 2^(scale-1);
        se = strel('disk',Nhood_rad,0);
        vf = mat2gray(sqrt(stdfilt(f,getnhood(se))));       
        im = cat(3,im,vf);
    end
    
    if (strcmp(which_features{feature_iter},'x'))

        %% x position feature
        xx = repmat((1:size(im_rgb,2)),[size(im_rgb,1),1]);
        xx = mat2gray(xx);
        im = cat(3,im,xx);
    end
    if (strcmp(which_features{feature_iter},'y'))

        %% y position feature
        yy = repmat((1:size(im_rgb,1))',[1 size(im_rgb,2)]);
        yy = mat2gray(yy);
        im = cat(3,im,yy);
    end
            
    if strcmp(which_features{feature_iter},'hue opp')
        % hue
        theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
        theta(white_index) = -pi/2;
        im_theta = reshape(theta,size(r));
        if opts.features.plot
            h1 = figure; imagesc(im_theta);
            cmap = colormap(hsv); colorbar('southoutside');
            axis equal; axis off; axis tight;set(gca,'FontSize',20);
            c1 = cmap(1:11,:); c2 = cmap(12:22,:); c3 = cmap(23:32,:);
            c4 = cmap(33:43,:); c5 = cmap(44:54,:); c6 = cmap(55:end,:);
            cmap_new = [c1;c2;c4;c5;c3;c6];
            h2 = figure; imagesc(im_theta);axis equal; axis off; axis tight;title('Hue');
            colormap(h1,cmap_new);colorbar('southoutside');set(gcf,'color','white');
            set(gca,'FontSize',20);
            close(h1);
        end
        %if ~opts.features.plot 
        %    close(h1); close(h2);
        %end
        
        im = im_theta;
    end
    
    if (strcmp(which_features{feature_iter},'saturation opp')) 
        % saturation
        sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
        im_sat = reshape(sat,size(r));
        if opts.features.plot
            figure; imagesc(im_sat); 
            colormap(jet); colorbar('southoutside'); %title('Saturation');
            axis equal; axis off; axis tight;set(gcf,'color','white');
            set(gca,'FontSize',20);
        end
        im = im_sat; 
    end
    
    if (strcmp(which_features{feature_iter},'brightness opp'))
        % brightness
        brightness = rotated_coordinates(1,:);
        im_brightness = reshape(brightness,size(r));
        if opts.features.plot
            figure; imagesc(im_brightness); 
            colormap(jet); colorbar('southoutside'); %title('Brightness');
            axis equal; axis off; axis tight;set(gcf,'color','white');
            set(gca,'FontSize',20);
        end
        im = im_brightness;
    end
    
    % downsample
    im = imresize(im,2^(-(scale-1)));
    
    %%
    if (opts.features.decorrelate) && ~ strcmp(which_features{feature_iter},'hue opp')...
            && ~ strcmp(which_features{feature_iter},'brightness opp')...
            && ~ strcmp(which_features{feature_iter},'saturation opp')
        im = mat2gray(pcaIm(im));
    end
    
%     % RESCALE THIS TO [0 255] 
%     if strcmp(opts.affinityFunction,'difference')
%     im = (im - min(im(:)))./(max(im(:)) - min(im(:)))*255;
%     im = uint8(im);
%     end
    %%
    f_maps{feature_iter} = im;
    end
end
