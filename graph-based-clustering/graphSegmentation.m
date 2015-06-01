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

function [segmented_image, E, E_oriented] = graphSegmentation(affinity_matrix,im_sizes,I,opts)
    % edge detection
    [E,E_oriented] = getE(affinity_matrix,im_sizes,I,opts);
    if opts.plot_results
        figure; imshow(1-mat2gray(E));
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
        E_ucm = contours2ucm_crisp_boundaries(mat2gray(E_oriented));
        %regionMerging; 
        segmented_image = ucm2colorsegs(E_ucm,I,thresh);toc;
        if opts.plot_results, figure;
        subplot(121); imshow(uint8(I)); subplot(122); 
        imshow(uint8(segmented_image)); end;
    else
        segmented_image = [];
    end
end
