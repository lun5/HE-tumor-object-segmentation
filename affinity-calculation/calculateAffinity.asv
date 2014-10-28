%% function [Ws,im_sizes] = getW(I,opts_affinity)
% builds affinity matrices Ws for image I
% 
% INPUTS
%  I                - NxMxC query image
%  opts_affinity    - parameter settings (see setEnvironment_affinity)
%
% OUTPUTS
%  Ws         - affinity matrices; Ws{i} is the affinity matrix for the image at scale i
%  im_sizes   - im_sizes{i} gives the dimensions of the image at scale i
%               (note: dimensions are num cols x num rows; this is the
%                opposite of matlab's default!)
%
% -------------------------------------------------------------------------
% HE segmentation Toolbox
% Luong Nguyen, 2014 [lun5@pitt.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------
% taken from 
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2014 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------
function [Ws,im_sizes] = calculateAffinity(I,opts_affinity)
    
    %%
    Ws = [];
    Ws_each_feature_set = [];
    im_sizes = [];
    num_scales = opts_affinity.num_scales;
    scale_offset = opts_affinity.scale_offset;
    
    %%
    for s=1:num_scales
        if (opts_affinity.display_progress), fprintf('\n\nProcessing scale %d:\n',s+scale_offset); end
        
        f_maps = [];
        for i=1:length(opts_affinity.features.which_features)
            f_maps{i} = getFeatures(double(I)/255,s+scale_offset,opts_affinity.features.which_features{i},opts_affinity);
        end
        
        %% NEED TO CHANGE HERE FOR DIFFERENT TYPE OF AFFINITY %%
        if strcmp(opts_affinity.affinityFunction,'PMI') 
        for feature_set_iter=1:length(f_maps)
            if (opts_affinity.display_progress), fprintf('\nProcessing feature type ''%s'':\n',opts_affinity.features.which_features{feature_set_iter}); end
        
            scale = 2^(-(s-1+scale_offset));
            
            f_maps_curr = f_maps{feature_set_iter};
            im_sizes{num_scales-s+1} = [size(f_maps_curr,2),size(f_maps_curr,1)];
            
            if ((s==1) || ~opts_affinity.only_learn_on_first_scale) % only learn models from first scale (and assume scale invariance)
                %% learn probability model
                if (opts_affinity.display_progress), fprintf('learning image model...'); tic; end
                p = learnP_A_B(f_maps_curr,opts_affinity);
                if (opts_affinity.display_progress), t = toc; fprintf('done: %1.2f sec\n', t); end
                %% learn w predictor
                if (opts_affinity.approximate_PMI)
                    if (opts_affinity.display_progress), fprintf('learning PMI predictor...'); tic; end
                    rf = learnPMIPredictor(f_maps_curr,p,opts_affinity);
                    if (opts_affinity.display_progress), t = toc; fprintf('done: %1.2f sec\n', t); end
                else
                    rf = [];
                end
            end
            
            %% build affinity matrix
            if (opts_affinity.display_progress), fprintf('building affinity matrix...'); tic; end
            if (strcmp(opts_affinity.model_type,'kde'))
                Ws_each_feature_set{num_scales-s+1}{feature_set_iter} = buildW_pmi(f_maps_curr,rf,p,opts_affinity);
            else
                error('unrecognized model type');
            end
            if (opts_affinity.display_progress), t = toc; fprintf('done: %1.2f sec\n', t); end
            
            %%
            if (feature_set_iter==1)
                Ws{num_scales-s+1} = Ws_each_feature_set{num_scales-s+1}{feature_set_iter};
            else
                Ws{num_scales-s+1} = Ws{num_scales-s+1}.*Ws_each_feature_set{num_scales-s+1}{feature_set_iter};
            end
        end
        end
    end    
end