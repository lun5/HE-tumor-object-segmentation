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
% Modified by Luong Nguyen 4/2015 for h&e hue channel
% 7/10/2015: calculate the join distribution at different sampling distance
% depending on parameter sigma in setEnvironment_affinity
function [Ws,Ws_each_feature_set, im_sizes] = getW(I,opts)
    
    %%
    Ws = [];
    Ws_each_feature_set = [];
    im_sizes = [];
    num_scales = opts.num_scales;
    scale_offset = opts.scale_offset;
    
    %%
    for s=1:num_scales
        if (opts.display_progress), fprintf('\n\nProcessing scale %d:\n',s+scale_offset); end
        
        
        [ f_maps, indx_white, indx_red] = getFeatures(double(I)./255,s+scale_offset,opts.features.which_features,opts);
        %num_pixels = size(I,1)*size(I,2);
        %% NEED TO CHANGE HERE FOR DIFFERENT TYPE OF AFFINITY %%
        
        for feature_set_iter=1:length(f_maps)
            if (opts.display_progress), fprintf('\nProcessing feature type ''%s'':\n',opts.features.which_features{feature_set_iter}); end
        
            scale = 2^(-(s-1+scale_offset));
            f_maps_curr = f_maps{feature_set_iter};
            nrows = size(f_maps_curr,1); ncols = size(f_maps_curr,2);
            num_pixels = nrows*ncols;
            im_sizes{num_scales-s+1} = [size(f_maps_curr,2),size(f_maps_curr,1)];
            which_feature = opts.features.which_features{feature_set_iter};
          if strcmp(opts.affinityFunction,'PMI') || strcmp(opts.affinityFunction,'PJoint')
            if (opts.display_progress), fprintf('\nProcessing affinity function ''%s'':\n',opts.affinityFunction); end
            if ((s==1) || ~opts.only_learn_on_first_scale) % only learn models from first scale (and assume scale invariance)
                %% learn probability model
                if (opts.display_progress), fprintf('learning image model...'); tic; end
                %% add different sampling distance using sigma
                sigma = opts.sig; sigma_values = [1];
                p = cell(length(sigma_values),1); 
                mixture_params = cell(length(sigma_values),1);
                rf = cell(length(sigma_values),1);
                for sigma_scale = 1:length(sigma_values)
                    opts.sig = sigma*sigma_values(sigma_scale);
                    if strcmp(which_feature,'hue opp')
                        Nsamples = opts.PMI_predictor.Nsamples_learning_PMI_predictor;
                        p{sigma_scale} = []; % we don't use kde to fit the joint distribution
                        %% I need to change here for different sigma
                        F = sampleF(f_maps_curr,Nsamples,opts);
                        %[ mu_hat_polar,~, kappa_hat,~, prior_probs] = moVM([cos(f_maps_curr(:)) sin(f_maps_curr(:))],3);
                        %init_params.prior_probs = prior_probs;
                        %init_params.theta_hat = mu_hat_polar;
                        %init_params.kappa_hat = kappa_hat;
                        Xcart = [cos(f_maps_curr(:)) sin(f_maps_curr(:))];
                        indx_purple_pink = (~indx_white) & (~indx_red);
                        Xcart_pp = Xcart(indx_purple_pink,:);
                        mu_white = 2.24; kappa_white = 30; %30 before 
                        [ mu_hat_polar_pp,~, kappa_hat_pp,~, prior_probs_pp] = moVM_fixWhite(Xcart_pp,2);
                        init_params.theta_hat = [mu_hat_polar_pp mu_white];
                        init_params.kappa_hat = [kappa_hat_pp kappa_white];
                        init_params.prior_probs = [prior_probs_pp(1:2)./num_pixels *sum(indx_purple_pink),...
                            sum(indx_white)/num_pixels, sum(indx_red)/num_pixels, prior_probs_pp(3)./num_pixels *sum(indx_purple_pink)];
                        if opts.model_half_space_only
                            [ params,~, prior_probs] = mixture_of_bivariate_VM(F, 6, init_params);
                        else
                            [ params,~, prior_probs] = mixture_of_bivariate_VM(F, 9, init_params);
                        end
                        mixture_params{sigma_scale} = struct('params',params,...
                            'prior_probs',prior_probs,'init_params',init_params);
                        %plotPMI_theta;
                    else %strcmp(which_feature,'luminance')
                        p{sigma_scale} = learnP_A_B(f_maps_curr,opts);
                        mixture_params{sigma_scale} = [];
                    end
                    %% learn w predictor
                    if (opts.approximate_PMI)
                        if (opts.display_progress), fprintf('learning PMI predictor...'); tic; end
                        rf{sigma_scale} = learnPMIPredictor(f_maps_curr,p{sigma_scale},mixture_params{sigma_scale}, which_feature, opts);
                    else
                        rf{sigma_scale} = [];
                    end
                    if (opts.display_progress), t = toc; fprintf('done: %1.2f sec\n', t); end
                    if (opts.display_progress), fprintf('building affinity matrix for feature %s at sigma = %1.1f ...',...
                            which_feature, opts.sig); tic; end
                    Ws_sampling_distance{num_scales-s+1}{feature_set_iter}{sigma_scale} = ...
                        buildW_pmi(f_maps_curr,rf{sigma_scale},p{sigma_scale},mixture_params{sigma_scale}, which_feature, opts);
                    if sigma_scale == 1
                        Ws_each_feature_set{num_scales-s+1}{feature_set_iter} = ...
                            Ws_sampling_distance{num_scales-s+1}{feature_set_iter}{sigma_scale};
                    else
                        Ws_each_feature_set{num_scales-s+1}{feature_set_iter} = ...
                            Ws_each_feature_set{num_scales-s+1}{feature_set_iter}.*Ws_sampling_distance{num_scales-s+1}{feature_set_iter}{sigma_scale};
                    end
                    if (opts.display_progress), t = toc; fprintf('done: %1.2f sec\n', t); end
                end
            end
          end
          
          %% build affinity matrix
          %if (opts.display_progress), fprintf('building affinity matrix...'); tic; end
          %Ws_each_feature_set{num_scales-s+1}{feature_set_iter} = ...
          %    buildW_pmi(f_maps_curr,rf,p,mixture_params, which_feature, opts);         
          
          %if (opts.display_progress), t = toc; fprintf('done: %1.2f sec\n', t); end
          
%             %% build affinity matrix
%             if (opts.display_progress), fprintf('building affinity matrix...'); tic; end
%             if (strcmp(opts.model_type,'kde'))
%                 Ws_each_feature_set{num_scales-s+1}{feature_set_iter} = buildW_pmi(f_maps_curr,rf,p,opts);
%             else 
%                 error('unrecognized model type');
%             end
%             if (opts.display_progress), t = toc; fprintf('done: %1.2f sec\n', t); end
            
%             elseif strcmp(opts.affinityFunction,'difference')
%                if (opts.display_progress), fprintf('\nProcessing affinity function ''%s'':\n',opts.affinityFunction); end
%                Ws_each_feature_set{num_scales-s+1}{feature_set_iter} = fast_buildW(f_maps_curr,opts.localPairs.rad, [], [], [],opts);
               %tElapsed = toc(tStart); fprintf('Affinity calculation takes: %1.2f sec\n',tElapsed);
         
             %%
            if (feature_set_iter==1)
                Ws{num_scales-s+1} = Ws_each_feature_set{num_scales-s+1}{feature_set_iter};
            else
                Ws{num_scales-s+1} = Ws{num_scales-s+1}.*Ws_each_feature_set{num_scales-s+1}{feature_set_iter};
            end
        end

    end  
    if opts.affinity.plot 
        %plotWs;
    end
end