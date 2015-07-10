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

function [opts] = setEnvironment_affinity
    %% scales                                                   used throughout code:
    opts.num_scales = 1;                               % how many image scales to measure affinity over
                                                                %  each subsequent scale is half the size of the one before (in both dimensions)
                                                                %  if opts.num_scales>1, then the method of Maire & Yu 2013 is used for globalization (going from affinty to boundaries);
                                                                %  otherwise regular spectral clustering is used
    opts.scale_offset = 0;                             % if opts.scale_offset==n then the first n image scales are skipped (first scales are highest resolution)
    
    
    %% features                                                 used in getFeatures.m:
    %opts.features.which_features = {'color','var'};            % which features to use:
    %opts.features.which_features = {'luminance'};
    %opts.features.which_features = {'hue opp', 'brightness opp', 'saturation opp'}; 
    opts.features.which_features = {'hue opp'};
    rotation_matrix = load(fullfile(pwd,'DanTrainingData','rotation_matrix_tp10-867-1.mat'),'rotation_matrix');
    opts.features.rotation_matrix = rotation_matrix.rotation_matrix;
    opts.features.decorrelate = 1;                              % decorrelate feature channels (done separately for each feature type in which_features)?
    opts.features.plot = false;
    %opts.features.rotation_matrix = rotation_matrix;
    %% Luong Nguyen 10/06/14 add opts.localPairs.rad,opts.localPairs.rad_inner used in 
    opts.localPairs.rad = 7;%5;
    opts.localPairs.rad_inner= [];
    opts.pyramid_ht = 1; % if we difference as a measure
    %% affinity function NEED TO INCLUDE THIS IN calculateAffinity 
    opts.affinityFunction = 'PMI';                           % PMI, differences, for now PMI
    %opts.affinityFunction = 'difference';     
    %% model and learning for PMI_{\rho}(A,B)                   used in learnP_A_B.m and buildW_pmi.m:
    opts.model_type = 'kde';                                    % what type of density estimate? (kde refers to kernel density estimation, which is the only method currently supported)
    opts.joint_exponent = 2; %1.25;                              % exponent \rho for PMI_{\rho} (Eqn. 2 in the paper)
    %opts.p_reg = 100;                                          % regularization added to numerator and demoninator of PMI calculation
    
    % kde options
    opts.kde.Nkernels = 10000;                                  % how many kernels for kde
    opts.kde.kdtree_tol = 0.001;                                % controls how exact is the kde evaluation (kde uses a kdtree to speed it up)
    opts.kde.learn_bw = true;                                   % adapt the bandwidth of the kde kernels to each test image?
    opts.kde.min_bw = 0.01; opts.kde.max_bw = 0.1;              % min and max bandwidths allowed when adapating bandwidth to test image
    
    % options for Eqn. 1 in paper
    opts.sig = 3;%0.5                                          % variance in pixels on Gaussian weighting function w(d) (see Eqn. 1 in paper)
    
    opts.model_half_space_only = false;                          % when true we model only half the joint {A,B} space and assume symmetry
    % speed up options
    opts.only_learn_on_first_scale = true;                      % setting this to true makes it so kde bandwidths and Affinity predictor are only 
                                                                %  learned on first scale (highest resolution) and assumed to be the same on lower 
                                                                %  resolution scales
                                                                %  (this often works well since image statistics are largely scale invariant)
    
                                                            
    %% approximate PMI with a random forest?                    used in learnPMIpredictor:                                
    opts.approximate_PMI = true;                                % approximate W with a random forest?
    opts.PMI_predictor.Nsamples_learning_PMI_predictor = 10000; % how many samples to learn approximation from
    opts.PMI_predictor.Ntrees = 4;                              % how many trees in the random forest
    
    %% display progress
    opts.display_progress = false;
    %% plot affinity matrix
    opts.affinity.plot = false;
end