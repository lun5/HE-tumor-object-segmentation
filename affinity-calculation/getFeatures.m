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

function [f_maps] = getFeatures(im_rgb,scale,which_feature,opts)
    
    im = [];
    
    %%
    if (strcmp(which_feature,'luminance'))
        im = mat2gray(mean(im_rgb,3));
    end
    if (strcmp(which_feature,'color'))
        %% color features
        if (size(im_rgb,3)==3)
            colorTransform = makecform('srgb2lab');
            cf = mat2gray(applycform(im_rgb./255,colorTransform));
            im = cat(3,im,cf(:,:,1:3));
        elseif (size(im_rgb,3)==1) % grayscale image
            im = im_rgb;
        else
            error('unhandled image format');
        end
        % NEED TO RESCALE THIS TO [0 255] 
    end
    
    if (strcmp(which_feature,'var'))    
        %% variance features
        f = pcaIm(im_rgb);
        
        Nhood_rad = 2^(scale-1);
        se = strel('disk',Nhood_rad,0);
        vf = mat2gray(sqrt(stdfilt(f,getnhood(se))));       
        im = cat(3,im,vf);
    end
    
    if (strcmp(which_feature,'x'))

        %% x position feature
        xx = repmat((1:size(im_rgb,2)),[size(im_rgb,1),1]);
        xx = mat2gray(xx);
        im = cat(3,im,xx);
    end
    if (strcmp(which_feature,'y'))

        %% y position feature
        yy = repmat((1:size(im_rgb,1))',[1 size(im_rgb,2)]);
        yy = mat2gray(yy);
        im = cat(3,im,yy);
    end
    
    if strcmp(which_feature,'hue opp') || (strcmp(which_feature,'saturation opp'))  || (strcmp(which_feature,'brightness opp'))
        getRotMat;
        opts.features.rotation_matrix = rotation_matrix;
        opts.features.decorrelate = 0;
        r = im_rgb(:,:,1)./255; g = im_rgb(:,:,2)./255; b = im_rgb(:,:,3)./255;
        rotated_coordinates = opts.features.rotation_matrix*double([r(:)'; g(:)'; b(:)']);
    end
    
    if strcmp(which_feature,'hue opp')
        % hue
        theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
        im_theta = reshape(theta,size(r));
        if opts.plot
            figure; imagesc(im_theta); 
            colormap(hsv); colorbar('southoutside'); title('Hue');
            axis equal; axis off; axis tight;
        end
        im = im_theta;
    end
    
    if (strcmp(which_feature,'saturation opp')) 
        % saturation
        sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
        im_sat = reshape(sat,size(r));
        if opts.plot
            figure; imagesc(im_sat); 
            colormap(jet); colorbar('southoutside'); title('Saturation');
            axis equal; axis off; axis tight;
        end
        im = im_sat; 
    end
    
    if (strcmp(which_feature,'brightness opp'))
        % brightness
        brightness = rotated_coordinates(1,:);
        im_brightness = reshape(brightness,size(r));
        if opts.plot
            figure; imagesc(im_brightness); 
            colormap(jet); colorbar('southoutside'); title('Brightness');
            axis equal; axis off; axis tight;
        end
        im = im_brightness;
    end
    
    %% downsample
    im = imresize(im,2^(-(scale-1)));
    
    %%
    if (opts.features.decorrelate)
        im = mat2gray(pcaIm(im));
    end
    
%     % RESCALE THIS TO [0 255] 
%     if strcmp(opts.affinityFunction,'difference')
%     im = (im - min(im(:)))./(max(im(:)) - min(im(:)))*255;
%     im = uint8(im);
%     end
    %%
    f_maps = im;
end