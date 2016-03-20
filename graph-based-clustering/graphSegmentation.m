% graphsegmentation
%% function [segmented_image, numComponents] = graphSegmentation(I, affinity_matrix, opts_clustering)
% given the image and affinity matrix
% 
% INPUTS
%  affinity_matrix   - affinity matrices; affinity_matrix{i} is the affinity matrix for the image at scale i
%  im_sizes   - im_sizes{i} gives the dimensions of the image at scale i
%               (note: dimensions are num cols x num rows; this is the
%                opposite of matlab's default!); Can I just get this from
%                affinity matrix?
%  I          - NxMxC query image
%  opts_clustering - parameter settings (see setEnvironment_clustering)
%
% OUTPUTS
%  E          - NxM boundary map
%  E_oriented - NxMxO boundary map split into boundaries energy at O orientations
% 
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2014 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

function [segmented_image_allfeatures,E_ucm_weighted, E_weighted, E_oriented] = graphSegmentation(affinity_matrix,aff_each_feature_set,im_sizes,I,opts)
    % edge detection
    num_features = length(aff_each_feature_set);
    E_ucm = []; segmented_image = [];
    for i = 1:num_features
        [E{i},E_oriented{i}] = getE({aff_each_feature_set{i}},im_sizes,I,opts); 
        E_oriented{i} = imresize(E_oriented{i},size(I(:,:,1)));
        E{i} = imresize(E{i},size(I(:,:,1)));
        if opts.plot_results
            figure; imshow(1-mat2gray(E{i}));
            %set(gca,'position',[0 0 1 1],'units','normalized')
            set(gcf,'color','white'); axis off; axis equal;axis tight;            
        end    
        %% Segment image
        % builds an Ultrametric Contour Map from the detected boundaries (E_oriented)
        % then segments image based on this map
        %
        % this part of the code is only supported on Mac and Linux
        if (~ispc) && opts.calculate_segments
            tic;thresh = 0.2;
            E_ucm{i} = contours2ucm_crisp_boundaries(mat2gray(E_oriented{i}));
            [segmented_image{i}, ~] = ucm2colorsegs(E_ucm{i},I,thresh);
            if opts.plot_results, figure; imshow(uint8(segmented_image{i}));end                          
        else
            segmented_image = [];
            E_ucm = [];
        end
    end
    E_ucm_weighted = [];
    E_weighted = [];
    segmented_image_allfeatures = [];
    
    if length(aff_each_feature_set) == 3
        max_val = max(max(max(max(E{1},E{2}),E{3})));        
        for i = 1:3
            E{i} = E{i}.*max_val./max(max(E{i}));
        end
        weights = [10 2 1]';
        E_allfeatures = cat(3,E{:});
        W = repmat(weights./sum(weights),1, size(E_allfeatures,1), size(E_allfeatures,2));
        W = permute(W,[2 3 1]);
        E_weighted = sum(E_allfeatures.*W,3);
        if (~ispc) && opts.calculate_segments
            E_ucm_allfeatures = cat(3,E_ucm{:});
            E_ucm_weighted = sum(E_ucm_allfeatures.*W,3);
            %[segmented_image_allfeatures,~] = ucm2colorsegs(E_ucm_weighted,I,thresh);
            if opts.plot_results, figure; imshow(uint8(segmented_image_allfeatures));end
        else
            E_ucm_weighted = [];
            segmented_image_allfeatures = [];
        end
    end
end
