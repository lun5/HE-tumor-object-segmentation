%% function [opts_affinity] = setEnvironment_affinity
% set parameter options for getting input images
% 
% INPUTS
%  fileType - svs or already tiled image 
%  tilesize - 
%  type - specifies which parameter set to use 
%          e.g., can take on values 'speedy' or 'accurate' 
%          feel free to define your custom types at end of this function
%  parallelize option
%  
% OUTPUTS
%  opts_input - selected parameters
% 
% -------------------------------------------------------------------------
% HE segmentation toolbox
% Luong Nguyen, 2014 [lun5@pitt.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

function [opts_affinity] = setEnvironment_affinity
    %% scales                                                   used throughout code:
    opts_affinity.num_scales = 3;                               % how many image scales to measure affinity over
                                                                %  each subsequent scale is half the size of the one before (in both dimensions)
                                                                %  if opts.num_scales>1, then the method of Maire & Yu 2013 is used for globalization (going from affinty to boundaries);
                                                                %  otherwise regular spectral clustering is used
    opts_affinity.scale_offset = 0;                             % if opts.scale_offset==n then the first n image scales are skipped (first scales are highest resolution)
    
    
    %% features                                                 used in getFeatures.m:
    opts_affinity.features.which_features = {'color','var'};            % which features to use:
    %opts_affinity.features.which_features = {'luminance'};
    %opts_affinity.features.which_features = {'hue opp'};%, 'brightness opp', 'saturation opp'}; 
    opts_affinity.features.decorrelate = 0;                              % decorrelate feature channels (done separately for each feature type in which_features)?
    opts_affinity.features.plot = true;
    %opts_affinity.features.rotation_matrix = rotation_matrix;
    %% Luong Nguyen 10/06/14 add opts.localPairs.rad,opts.localPairs.rad_inner used in 
    opts_affinity.localPairs.rad = 5;
    opts_affinity.localPairs.rad_inner= [];
    
    %% affinity function NEED TO INCLUDE THIS IN calculateAffinity 
    opts_affinity.affinityFunction = 'PMI';                           % PMI, differences, for now PMI
%   opts_affinity.affinityFunction = 'difference';     
    %% model and learning for PMI_{\rho}(A,B)                   used in learnP_A_B.m and buildW_pmi.m:
    opts_affinity.model_type = 'kde';                                    % what type of density estimate? (kde refers to kernel density estimation, which is the only method currently supported)
    opts_affinity.joint_exponent = 1.25;                                 % exponent \rho for PMI_{\rho} (Eqn. 2 in the paper)
    opts_affinity.p_reg = 100;                                           % regularization added to numerator and demoninator of PMI calculation
    
    % kde options
    opts_affinity.kde.Nkernels = 1000;                                  % how many kernels for kde
    opts_affinity.kde.kdtree_tol = 0.001;                                % controls how exact is the kde evaluation (kde uses a kdtree to speed it up)
    opts_affinity.kde.learn_bw = false;                                   % adapt the bandwidth of the kde kernels to each test image?
    opts_affinity.kde.min_bw = 0.01; opts_affinity.kde.max_bw = 0.1;              % min and max bandwidths allowed when adapating bandwidth to test image
    
    % options for Eqn. 1 in paper
    opts_affinity.sig = 0.25;                                            % variance in pixels on Gaussian weighting function w(d) (see Eqn. 1 in paper)
    
    % speed up options
    opts_affinity.only_learn_on_first_scale = true;            % setting this to true makes it so kde bandwidths and Affinity predictor are only 
                                                                %  learned on first scale (highest resolution) and assumed to be the same on lower 
                                                                %  resolution scales
                                                                %  (this often works well since image statistics are largely scale invariant)
    
                                                            
    %% approximate PMI with a random forest?                    used in learnPMIpredictor:                                
    opts_affinity.approximate_PMI = true;                                % approximate W with a random forest?
    opts_affinity.PMI_predictor.Nsamples_learning_PMI_predictor = 1000; % how many samples to learn approximation from
    opts_affinity.PMI_predictor.Ntrees = 4;                              % how many trees in the random forest
    
    %% display progress
    opts_affinity.display_progress = true;
    %% plot affinity matrix
    opts_affinity.plot = true;
end